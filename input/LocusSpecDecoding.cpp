//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "input/LocusSpecDecoding.hh"

#include <memory>

#include <boost/optional.hpp>

#include "input/GraphBlueprint.hh"
#include "input/RegionGraph.hh"

using boost::optional;
using graphtools::Graph;
using graphtools::NodeId;
using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

static vector<GenomicRegion> addFlankingRegions(int kExtensionLength, const vector<GenomicRegion>& referenceRegions)
{
    const GenomicRegion& firstRegion = referenceRegions.front();
    const int64_t leftFlankStart = firstRegion.start() - kExtensionLength;
    GenomicRegion leftFlank(firstRegion.contigIndex(), leftFlankStart, firstRegion.start());

    const GenomicRegion& lastRegion = referenceRegions.back();
    const int64_t rightFlankEnd = lastRegion.end() + kExtensionLength;
    GenomicRegion rightFlank(lastRegion.contigIndex(), lastRegion.end(), rightFlankEnd);

    auto regions = referenceRegions;
    regions.insert(regions.begin(), std::move(leftFlank));
    regions.push_back(std::move(rightFlank));

    return regions;
}

static string extendLocusStructure(
    const Reference& reference, const vector<GenomicRegion>& referenceRegions,
    const string& flanklessLocusStructure)
{

    const auto& leftFlank = referenceRegions.front();
    const auto& rightFlank = referenceRegions.back();
    return reference.getSequence(leftFlank) + flanklessLocusStructure
        + reference.getSequence(rightFlank);
}

static vector<GenomicRegion>
addReferenceRegionsForInterruptions(const GraphBlueprint& blueprint, const vector<GenomicRegion>& referenceRegions)
{
    int regionIndex = 0;
    vector<GenomicRegion> completedReferenceRegions;

    for (const auto& feature : blueprint)
    {
        if (feature.type == GraphBlueprintFeatureType::kInterruption)
        {
            assert(regionIndex != 0 && regionIndex < static_cast<int>(referenceRegions.size()));
            const auto& leftRegion = referenceRegions[regionIndex];
            const auto& rightRegion = referenceRegions[regionIndex + 1];
            completedReferenceRegions.emplace_back(leftRegion.contigIndex(), leftRegion.end(), rightRegion.start());
        }
        else
        {
            completedReferenceRegions.push_back(referenceRegions[regionIndex]);
            ++regionIndex;
        }
    }

    assert(blueprint.size() == completedReferenceRegions.size());
    return completedReferenceRegions;
}

static vector<string>
combine(const std::string& prefix, const ReferenceContigInfo& contigInfo, const vector<GenomicRegion>& regions)
{
    vector<string> combinedStrings;
    for (const auto& region : regions)
    {
        combinedStrings.push_back(prefix + "_" + encode(contigInfo, region));
    }

    return combinedStrings;
}

static GenomicRegion mergeRegions(const vector<GenomicRegion>& regions)
{
    const int kMaxMergeDistance = 500;
    vector<GenomicRegion> mergedReferenceRegions = merge(regions, kMaxMergeDistance);
    if (mergedReferenceRegions.size() != 1)
    {
        std::stringstream out;
        for (const GenomicRegion& region : regions)
        {
            out << region << " ";
        }
        throw std::runtime_error(
            "Expected reference regions to be closer than " + to_string(kMaxMergeDistance)
            + " from one another: " + out.str());
    }

    return mergedReferenceRegions.front();
}

AlleleCount determineExpectedAlleleCount(Sex sex, const string& chrom)
{
    const bool isFemaleChromY = sex == Sex::kFemale && (chrom == "chrY" || chrom == "Y");
    if (isFemaleChromY)
    {
        return AlleleCount::kZero;
    }

    const bool isSexChrom = chrom == "chrX" || chrom == "X" || chrom == "chrY" || chrom == "Y";
    if (sex == Sex::kMale && isSexChrom)
    {
        return AlleleCount::kOne;
    }

    return AlleleCount::kTwo;
}

static NodeToRegionAssociation associateNodesWithReferenceRegions(
    const GraphBlueprint& blueprint, const Graph& graph, const vector<GenomicRegion>& referenceRegions)
{
    assert(blueprint.size() == referenceRegions.size());

    NodeToRegionAssociation referenceRegionsOfGraphNodes;

    for (int featureIndex = 0; featureIndex != static_cast<int>(blueprint.size()); ++featureIndex)
    {
        const auto& feature = blueprint[featureIndex];
        const auto& referenceRegion = referenceRegions[featureIndex];

        for (const auto& nodeId : feature.nodeIds)
        {
            const int nodeLength = graph.nodeSeq(nodeId).length();
            GenomicRegion referenceRegionForNode(
                referenceRegion.contigIndex(), referenceRegion.start(), referenceRegion.start() + nodeLength);
            referenceRegionsOfGraphNodes.emplace(std::make_pair(nodeId, std::move(referenceRegionForNode)));
        }
    }

    return referenceRegionsOfGraphNodes;
}

static VariantType determineVariantType(GraphBlueprintFeatureType featureType)
{
    switch (featureType)
    {
    case GraphBlueprintFeatureType::kInsertionOrDeletion:
    case GraphBlueprintFeatureType::kSwap:
        return VariantType::kSmallVariant;
    case GraphBlueprintFeatureType::kSkippableRepeat:
    case GraphBlueprintFeatureType::kUnskippableRepeat:
        return VariantType::kRepeat;
    default:
        std::ostringstream encoding;
        encoding << featureType;
        throw std::logic_error("Feature of type " + encoding.str() + " does not define a variant");
    }
}

static VariantSubtype determineVariantSubtype(
    GraphBlueprintFeatureType featureType, VariantTypeFromUser userDescription, const GenomicRegion referenceRegion)
{
    if (featureType == GraphBlueprintFeatureType::kInsertionOrDeletion)
    {
        if (referenceRegion.length() == 0)
        {
            return VariantSubtype::kInsertion;
        }
        else
        {
            return VariantSubtype::kDeletion;
        }
    }
    else if (featureType == GraphBlueprintFeatureType::kSwap)
    {
        if (userDescription == VariantTypeFromUser::kSMN)
        {
            return VariantSubtype::kSMN;
        }
        else
        {
            return VariantSubtype::kSwap;
        }
    }
    else if (userDescription == VariantTypeFromUser::kCommonRepeat)
    {
        return VariantSubtype::kCommonRepeat;
    }
    else if (userDescription == VariantTypeFromUser::kRareRepeat)
    {
        return VariantSubtype::kRareRepeat;
    }
    else
    {
        std::ostringstream message;
        message << featureType;
        throw std::logic_error("Feature " + message.str() + " does not correspond to variant");
    }
}

static optional<NodeId> determineReferenceNode(
    const GraphBlueprintFeature& feature,
    const Reference& reference,
    const GenomicRegion& referenceRegion)
{

    if (feature.type == GraphBlueprintFeatureType::kSkippableRepeat
        || feature.type == GraphBlueprintFeatureType::kUnskippableRepeat)
    {
        return feature.nodeIds.front();
    }

    const auto& contigName = reference.contigInfo().getContigName(referenceRegion.contigIndex());
    const string refSequence = reference.getSequence(contigName, referenceRegion.start(), referenceRegion.end());

    optional<NodeId> optionalReferenceNode;
    for (int index = 0; index != static_cast<int>(feature.nodeIds.size()); ++index)
    {
        if (refSequence == feature.sequences[index])
        {
            optionalReferenceNode = feature.nodeIds[index];
            break;
        }
    }

    return optionalReferenceNode;
}

LocusSpecification decodeLocusSpecification(
    const LocusDescriptionFromUser& userDescription,
    Sex sampleSex,
    const Reference& reference,
    const HeuristicParameters& heuristicParams)
{
    assertValidity(userDescription);

    const int kExtensionLength = heuristicParams.regionExtensionLength();
    auto referenceRegionsWithFlanks = addFlankingRegions(kExtensionLength, userDescription.referenceRegions);
    auto completeLocusStructure = extendLocusStructure(
        reference, referenceRegionsWithFlanks, userDescription.locusStructure);

    GraphBlueprint blueprint = decodeFeaturesFromRegex(completeLocusStructure);
    graphtools::Graph locusGraph = makeRegionGraph(blueprint, userDescription.locusId);
    auto completeReferenceRegions = addReferenceRegionsForInterruptions(blueprint, referenceRegionsWithFlanks);

    vector<string> variantIds
        = combine(userDescription.locusId, reference.contigInfo(), userDescription.referenceRegions);

    GenomicRegion referenceRegionForEntireLocus = mergeRegions(userDescription.referenceRegions);

    vector<GenomicRegion> targetReadExtractionRegions;
    for (const GenomicRegion& region : userDescription.targetRegions)
    {
        targetReadExtractionRegions.push_back(region.extend(kExtensionLength));
    }
    if (targetReadExtractionRegions.empty())
    {
        targetReadExtractionRegions.push_back(referenceRegionForEntireLocus.extend(kExtensionLength));
    }

    const auto& contigName
        = reference.contigInfo().getContigName(referenceRegionForEntireLocus.contigIndex());
    AlleleCount expectedAlleleCount = determineExpectedAlleleCount(sampleSex, contigName);

    NodeToRegionAssociation referenceRegionsOfGraphNodes
        = associateNodesWithReferenceRegions(blueprint, locusGraph, completeReferenceRegions);

    LocusSpecification regionSpec(
        userDescription.locusId, targetReadExtractionRegions, expectedAlleleCount, locusGraph,
        referenceRegionsOfGraphNodes);
    regionSpec.setOfftargetReadExtractionRegions(userDescription.offtargetRegions);

    int variantIndex = 0;
    for (const auto& feature : blueprint)
    {
        if (doesFeatureDefineVariant(feature.type))
        {
            const GenomicRegion& referenceRegion = userDescription.referenceRegions.at(variantIndex);

            VariantTypeFromUser variantDescription = userDescription.variantTypesFromUser.at(variantIndex);
            VariantType variantType = determineVariantType(feature.type);
            VariantSubtype variantSubtype = determineVariantSubtype(feature.type, variantDescription, referenceRegion);

            optional<NodeId> optionalReferenceNode
                = determineReferenceNode(feature, reference, referenceRegion);

            VariantClassification classification(variantType, variantSubtype);

            regionSpec.addVariantSpecification(
                variantIds[variantIndex], classification, referenceRegion, feature.nodeIds, optionalReferenceNode);

            ++variantIndex;
        }
    }

    return regionSpec;
}

void assertValidity(const LocusDescriptionFromUser& userDescription)
{
    const GraphBlueprint blueprint = decodeFeaturesFromRegex(userDescription.locusStructure);
    int numVariants = 0;
    for (const GraphBlueprintFeature& feature : blueprint)
    {
        if (doesFeatureDefineVariant(feature.type))
        {
            ++numVariants;
        }
    }

    if (numVariants == 0)
    {
        throw std::runtime_error(
            "Locus " + userDescription.locusId + " must encode at least one variant " + userDescription.locusStructure);
    }

    if (numVariants != static_cast<int>(userDescription.referenceRegions.size()))
    {
        throw std::runtime_error(
            "Locus " + userDescription.locusId + " must specify reference regions for " + to_string(numVariants)
            + " variants");
    }

    if (numVariants != static_cast<int>(userDescription.variantTypesFromUser.size()))
    {
        throw std::runtime_error(
            "Locus " + userDescription.locusId + " must specify variant types for " + to_string(numVariants)
            + " variants");
    }
}

}
