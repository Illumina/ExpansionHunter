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

#include "input/LocusSpecDecoding.hh"

#include <algorithm>
#include <memory>
#include <stdexcept>

#include <boost/optional.hpp>

#include "common/WorkflowContext.hh"
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

static CopyNumberBySex determineCopyNumber(const string& contig)
{
    if (contig == "chrY" || contig == "Y")
    {
        return CopyNumberBySex::kZeroInFemaleOneInMale;
    }

    if (contig == "chrX" || contig == "X")
    {
        return CopyNumberBySex::kTwoInFemaleOneInMale;
    }

    return CopyNumberBySex::kTwoInFemaleTwoInMale;
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

static VariantSubtype determineCnvVariantSubtype(VariantSubtypeFromUser variantSubtypeFromUser)
{
    if (variantSubtypeFromUser == VariantSubtypeFromUser::kTarget)
    {
        return VariantSubtype::kTarget;
    }
    if (variantSubtypeFromUser == VariantSubtypeFromUser::kBaseline)
    {
        return VariantSubtype::kBaseline;
    }
    else
    {
        throw std::logic_error("Illegal CNV variant subtype");
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

static vector<GenomicRegion> computeStatsRegions(const GenomicRegion& locusLocation, int flankLength)
{
    const int64_t leftFlankStart = locusLocation.start() - flankLength;
    const int64_t leftFlankEnd = locusLocation.start();
    GenomicRegion leftFlank(locusLocation.contigIndex(), leftFlankStart, leftFlankEnd);

    const int64_t rightFlankStart = locusLocation.end();
    const int64_t rightFlankEnd = locusLocation.end() + flankLength;
    GenomicRegion rightFlank(locusLocation.contigIndex(), rightFlankStart, rightFlankEnd);

    return { leftFlank, rightFlank };
}

GraphLocusSpec
decodeGraphLocusSpecification(const LocusDescriptionFromUser& userDescription, const Reference& reference)
{
    try
    {
        assertValidity(userDescription);

        WorkflowContext context;

        vector<GenomicRegion> variantLocations;
        for (const auto& variant : userDescription.variantDescriptionFromUsers)
        {
            variantLocations.push_back(variant.variantLocation);
        }
        const int kExtensionLength = context.heuristics().regionExtensionLength();
        auto referenceRegionsWithFlanks = addFlankingRegions(kExtensionLength, variantLocations);
        assert(userDescription.locusStructure);
        auto locusStructure = userDescription.locusStructure;
        auto completeLocusStructure = extendLocusStructure(reference, referenceRegionsWithFlanks, *locusStructure);

        const auto& locusId = userDescription.locusId;
        GraphBlueprint blueprint = decodeFeaturesFromRegex(completeLocusStructure);
        graphtools::Graph locusGraph = makeRegionGraph(blueprint, locusId);
        auto completeReferenceRegions = addReferenceRegionsForInterruptions(blueprint, referenceRegionsWithFlanks);

        vector<GenomicRegion> targetReadExtractionRegions;
        for (const GenomicRegion& region : userDescription.targetRegions)
        {
            targetReadExtractionRegions.push_back(region.extend(kExtensionLength));
        }
        if (targetReadExtractionRegions.empty())
        {

            targetReadExtractionRegions.push_back(userDescription.locusLocation.extend(kExtensionLength));
        }

        const auto& contigName = reference.contigInfo().getContigName(userDescription.locusLocation.contigIndex());
        auto copyNumber = determineCopyNumber(contigName);

        NodeToRegionAssociation referenceRegionsOfGraphNodes
            = associateNodesWithReferenceRegions(blueprint, locusGraph, completeReferenceRegions);

        GenotyperParameters parameters;
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

        GraphLocusReferenceRegions referenceRegions;
        referenceRegions.offtargetRegionsWithReads = userDescription.offtargetRegions;
        referenceRegions.regionsWithReads = targetReadExtractionRegions;
        referenceRegions.statsRegions = computeStatsRegions(userDescription.locusLocation, kExtensionLength);

        ReferenceGraph referenceGraph(locusGraph, referenceRegionsOfGraphNodes);

        GraphLocusSpec locusSpec(locusId, copyNumber, referenceRegions, referenceGraph, parameters);

        int variantIndex = 0;
        for (const auto& feature : blueprint)
        {
            if (doesFeatureDefineVariant(feature.type))
            {
                VariantDescriptionFromUser variant = userDescription.variantDescriptionFromUsers.at(variantIndex);
                const GenomicRegion& referenceRegion = variant.variantLocation;

                VariantTypeFromUser variantTypeFromDescription = variant.variantType;
                const string& variantId = variant.variantId;
                VariantType variantType = determineVariantType(feature.type);
                VariantSubtype variantSubtype
                    = determineVariantSubtype(feature.type, variantTypeFromDescription, referenceRegion);

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

CnvLocusSpec decodeCnvLocusSpecification(const LocusDescriptionFromUser& userDescription, const Reference& reference)
{
    try
    {
        const auto& contigName = reference.contigInfo().getContigName(userDescription.locusLocation.contigIndex());
        auto copyNumber = determineCopyNumber(contigName);
        CnvType cnvType = CnvType::kNonoverlapping;

        GenotyperParameters parameters;
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

        for (const auto& variant : userDescription.variantDescriptionFromUsers)
        {
            if (variant.variantSubtype == VariantSubtypeFromUser::kBaseline && !(*variant.expectedNormalCN))
            {
                cnvType = CnvType::kOverlapping;
            }
        }

        CnvLocusSpec locusSpec(userDescription.locusId, cnvType, copyNumber, parameters);

        for (const auto& variant : userDescription.variantDescriptionFromUsers)
        {
            const GenomicRegion& referenceRegion = variant.variantLocation;
            const string& variantId = variant.variantId;
            VariantType variantType = VariantType::kCNV;
            VariantSubtype variantSubtype = determineCnvVariantSubtype(*variant.variantSubtype);

            VariantClassification classification(variantType, variantSubtype);

            CnvGenotyperParameters variantParameters;
            variantParameters.regionGC = *variant.regionGC;
            variantParameters.mappingQualityThreshold = *variant.mappingQualityThreshold;
            variantParameters.maxCopyNumber = *variant.maxCopyNumber;
            variantParameters.depthScaleFactor = *variant.depthScaleFactor;
            variantParameters.standardDeviationOfCN2 = *variant.standardDevidationOfCN2;
            variantParameters.meanDepthValues = *variant.meanDepthValues;
            variantParameters.priorCopyNumberFrequency = *variant.priorCopyNumberFrequency;
            variantParameters.expectedNormal = *variant.expectedNormalCN;

            locusSpec.addVariantSpecification(variantId, classification, referenceRegion, variantParameters);
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
    assert(userDescription.locusStructure);
    auto locusStructure = userDescription.locusStructure;
    const GraphBlueprint blueprint = decodeFeaturesFromRegex(*locusStructure);
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
            "Locus " + userDescription.locusId + " must encode at least one variant " + *locusStructure);
    }

    if (numVariants != static_cast<int>(userDescription.variantDescriptionFromUsers.size()))
    {
        throw std::runtime_error(
            "Locus " + userDescription.locusId + " must specify variant information for " + to_string(numVariants)
            + " variants");
    }
}
}
