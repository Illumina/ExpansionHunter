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

#include "locus_spec/Decoding.hh"

#include "common/WorkflowContext.hh"
#include "locus_spec/GraphBlueprint.hh"
#include "locus_spec/RegionGraph.hh"

namespace ehunter
{

using std::string;
using std::to_string;
using std::vector;

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

/*
 * struct
{
    // Regions in the reference where we expect relevant reads to align
    std::vector<GenomicRegion> regionsWithReads;
    // Regions where additional relevant reads might be found that require filtering or special considerations
    std::vector<GenomicRegion> offtargetRegionsWithReads;
    std::vector<GenomicRegion> statsRegions;
};
 */

static vector<GenomicRegion> decodeRegions(const Reference& reference, const vector<string>& encodings)
{
    const ReferenceContigInfo& contigInfo = reference.contigInfo();

    vector<GenomicRegion> regions;
    regions.reserve(encodings.size());
    for (const auto& encoding : encodings)
    {
        regions.push_back(decode(contigInfo, encoding));
    }
    return regions;
}

static vector<GenomicRegion> computeFlanks(const vector<GenomicRegion>& regions, int extensionLength)
{
    vector<GenomicRegion> flanks;

    for (const auto& region : regions)
    {
        int64_t leftFlankStart = region.start() - extensionLength;
        int64_t leftFlankEnd = region.start();
        flanks.emplace_back(region.contigIndex(), leftFlankStart, leftFlankEnd);

        int64_t rightFlankStart = region.end();
        int64_t rightFlankEnd = region.end() + extensionLength;
        flanks.emplace_back(region.contigIndex(), rightFlankStart, rightFlankEnd);
    }

    return flanks;
}

static string addFlanks(
    const Reference& reference, const string& locusStructure, const GenomicRegion& locusLocation, int extensionLength)
{
    auto flanks = computeFlanks({ locusLocation }, extensionLength);

    string leftFlank = reference.getSequence(flanks.front());
    string rightFlank = reference.getSequence(flanks.back());

    const size_t maxNsAllowedInFlanks = 5;
    const size_t numNsInLeftFlank = std::count(leftFlank.begin(), leftFlank.end(), 'N');
    const size_t numNsInRightFlank = std::count(rightFlank.begin(), rightFlank.end(), 'N');

    if (numNsInLeftFlank + numNsInRightFlank > maxNsAllowedInFlanks)
    {
        const string message = "Flanks can contain at most " + to_string(maxNsAllowedInFlanks)
            + " characters N but found " + to_string(numNsInLeftFlank + numNsInRightFlank) + " Ns";
        throw std::runtime_error(message);
    }

    return leftFlank + locusStructure + rightFlank;
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

/*
static NodeLocations associateNodesWithReferenceRegions(
    const GraphBlueprint& blueprint, const Graph& graph, const vector<GenomicRegion>& referenceRegions)
{
    assert(blueprint.size() == referenceRegions.size());

    NodeLocations referenceRegionsOfGraphNodes;

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
} */

GraphLocusSpec decode(const Reference& reference, const GraphLocusEncoding& locusEncoding)
{
    // assertValidity(userDescription);

    // WorkflowContext context;
    // int flankLength = context.heuristics().regionExtensionLength();

    // auto locusStructureWithFlanks = addFlanks(reference, encoding.locusStructure, locusLocation, extensionLength);
    // GraphBlueprint blueprint = decodeFeaturesFromRegex(locusStructureWithFlanks);

    // GraphBlueprint blueprint = decode(locusEncoding);

    // const ReferenceContigInfo& contigInfo = reference.contigInfo();
    // auto locusLocation = decode(contigInfo, encoding.locusLocation);
    // auto contigCopyNumber = determineCopyNumber(contigInfo.getContigName(locusLocation.contigIndex()));

    /*
GraphLocusReferenceRegions referenceRegions;
referenceRegions.regionsWithReads = decodeRegions(reference, encoding.targetRegions);
referenceRegions.offtargetRegionsWithReads = decodeRegions(reference, encoding.offtargetRegions);
referenceRegions.statsRegions = computeFlanks(referenceRegions.regionsWithReads, extensionLength); */

    // graphtools::Graph locusGraph = makeRegionGraph(blueprint, encoding.locusId);

    // NodeLocations nodeLocations = getReferenceLocations(blueprint, locusLocation, encoding.variants);

    // auto completeReferenceRegions = addReferenceRegionsForInterruptions(blueprint, referenceRegionsWithFlanks);
    // NodeLocations referenceRegionsOfGraphNodes
    //    = associateNodesWithReferenceRegions(blueprint, locusGraph, completeReferenceRegions);
    // ReferenceGraph referenceGraph(locusGraph, referenceRegionsOfGraphNodes);

    GraphLocusSpec locusSpec(encoding.locusId, contigCopyNumber, referenceRegions, referenceGraph, genotyperParams);

    //

    /*
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

    NodeLocations referenceRegionsOfGraphNodes
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
    } */

    return locusSpec;
}

}