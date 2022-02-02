//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//

#include "io/LocusSpecDecoding.hh"

#include <algorithm>
#include <memory>
#include <stdexcept>

#include <boost/optional.hpp>

#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

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
    const Reference& reference, const vector<GenomicRegion>& referenceRegions, const string& flanklessLocusStructure)
{

    const auto& leftFlankRegion = referenceRegions.front();
    string leftFlank = reference.getSequence(leftFlankRegion);

    const auto& rightFlankRegion = referenceRegions.back();
    string rightFlank = reference.getSequence(rightFlankRegion);

    const size_t maxNsAllowedInFlanks = 5;
    const size_t numNsInLeftFlank = std::count(leftFlank.begin(), leftFlank.end(), 'N');
    const size_t numNsInRightFlank = std::count(rightFlank.begin(), rightFlank.end(), 'N');

    if (numNsInLeftFlank + numNsInRightFlank > maxNsAllowedInFlanks)
    {
        const string message = "Flanks can contain at most " + to_string(maxNsAllowedInFlanks)
            + " characters N but found " + to_string(numNsInLeftFlank + numNsInRightFlank) + " Ns";
        throw std::runtime_error(message);
    }

    return leftFlank + flanklessLocusStructure + rightFlank;
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

static ChromType determineChromosomeType(const string& chrom)
{
    if (chrom == "chrY" || chrom == "Y")
    {
        return ChromType::kY;
    }

    if (chrom == "chrX" || chrom == "X")
    {
        return ChromType::kX;
    }

    return ChromType::kAutosome;
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
    const GraphBlueprintFeature& feature, const Reference& reference, const GenomicRegion& referenceRegion)
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
    const LocusDescriptionFromUser& userDescription, const Reference& reference,
    const HeuristicParameters& heuristicParams)
{
    try
    {
        assertValidity(userDescription);

        const int kExtensionLength = heuristicParams.regionExtensionLength();
        auto referenceRegionsWithFlanks = addFlankingRegions(kExtensionLength, userDescription.referenceRegions);
        auto completeLocusStructure
            = extendLocusStructure(reference, referenceRegionsWithFlanks, userDescription.locusStructure);

        GraphBlueprint blueprint = decodeFeaturesFromRegex(completeLocusStructure);
        graphtools::Graph locusGraph = makeRegionGraph(blueprint, userDescription.locusId);
        auto completeReferenceRegions = addReferenceRegionsForInterruptions(blueprint, referenceRegionsWithFlanks);

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

        const auto& contigName = reference.contigInfo().getContigName(referenceRegionForEntireLocus.contigIndex());
        ChromType chromType = determineChromosomeType(contigName);

        NodeToRegionAssociation referenceRegionsOfGraphNodes
            = associateNodesWithReferenceRegions(blueprint, locusGraph, completeReferenceRegions);

        GenotyperParameters parameters(heuristicParams.minLocusCoverage());
        if (userDescription.errorRate)
        {
            parameters.errorRate = *userDescription.errorRate;
        }
        if (userDescription.likelihoodRatioThreshold)
        {
            parameters.likelihoodRatioThreshold = *userDescription.likelihoodRatioThreshold;
        }
        if (userDescription.minLocusCoverage)
        {
            parameters.minLocusCoverage = *userDescription.minLocusCoverage;
        }

        LocusSpecification locusSpec(
            userDescription.locusId, chromType, std::move(targetReadExtractionRegions), std::move(locusGraph),
            std::move(referenceRegionsOfGraphNodes), std::move(parameters), userDescription.useRFC1MotifAnalysis);
        locusSpec.setOfftargetReadExtractionRegions(userDescription.offtargetRegions);

        int variantIndex = 0;
        for (const auto& feature : blueprint)
        {
            if (doesFeatureDefineVariant(feature.type))
            {
                const GenomicRegion& referenceRegion = userDescription.referenceRegions.at(variantIndex);

                VariantTypeFromUser variantDescription = userDescription.variantTypesFromUser.at(variantIndex);
                const string& variantId = userDescription.variantIds[variantIndex];
                VariantType variantType = determineVariantType(feature.type);
                VariantSubtype variantSubtype
                    = determineVariantSubtype(feature.type, variantDescription, referenceRegion);

                optional<NodeId> optionalReferenceNode = determineReferenceNode(feature, reference, referenceRegion);

                VariantClassification classification(variantType, variantSubtype);

                locusSpec.addVariantSpecification(
                    variantId, classification, referenceRegion, feature.nodeIds, optionalReferenceNode);

                ++variantIndex;
            }
        }
        return locusSpec;
    }
    catch (const std::exception& e)
    {
        throw std::runtime_error("Error loading locus " + userDescription.locusId + ": " + e.what());
    }
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

    if (userDescription.useRFC1MotifAnalysis)
    {
        if ((numVariants != 1) or (userDescription.variantTypesFromUser[0] != VariantTypeFromUser::kCommonRepeat))
        {
            throw std::runtime_error(
                "Locus " + userDescription.locusId
                + " has option 'useRFC1MotifAnalysis' enabled, which requires that"
                  " exactly one variant of type 'Repeat' is defined.");
        }
    }
}

}
