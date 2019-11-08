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

#include "locus_spec/GraphLocusDecoding.hh"

#include <stdexcept>

#include <boost/optional.hpp>

#include "graphcore/Graph.hh"

#include "locus_spec/GraphBlueprint.hh"
#include "locus_spec/RegionGraph.hh"

using boost::optional;
using graphtools::Graph;
using graphtools::NodeId;
using std::runtime_error;
using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

static GenomicRegion getLocusLocation(const GraphLocusDecoding& locusEncoding)
{
    vector<GenomicRegion> variantLocations;
    for (const auto& variant : locusEncoding.variants)
    {
        variantLocations.push_back(variant.location);
    }

    const int kMaxMergeDistance = 500;
    vector<GenomicRegion> mergedLocations = merge(variantLocations, kMaxMergeDistance);
    if (mergedLocations.size() != 1)
    {
        throw runtime_error("Distance between variants must not exceed " + to_string(kMaxMergeDistance) + "bp");
    }

    return mergedLocations.front();
}

static CopyNumberBySex getCopyNumber(const string& contig)
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

static AnalysisRegions getAnalysisRegions(const GraphLocusDecoding& encoding, const GenomicRegion& locusLocation)
{
    AnalysisRegions regions;
    for (auto region : encoding.targetRegions)
    {
        regions.regionsWithReads.push_back(region.extend(encoding.flankLength));
    }
    for (auto variant : encoding.variants)
    {
        regions.regionsWithReads.push_back(variant.location.extend(encoding.flankLength));
    }

    regions.offtargetRegionsWithReads = encoding.offtargetRegions;

    // TODO: Handle case when user explicitly defines stats regions
    // TODO: Decide if flank length is always appropriate for stats regions
    regions.statsRegions = computeFlanks({ locusLocation }, encoding.flankLength);
    return regions;
}

static string
addFlanks(const Reference& reference, const string& structure, const GenomicRegion& location, int flankLength)
{
    auto flanks = computeFlanks({ location }, flankLength);
    string leftFlank = reference.getSequence(flanks.front());
    string rightFlank = reference.getSequence(flanks.back());

    const size_t maxNsAllowed = 5;
    const size_t numNsInLeftFlank = std::count(leftFlank.begin(), leftFlank.end(), 'N');
    const size_t numNsInRightFlank = std::count(rightFlank.begin(), rightFlank.end(), 'N');
    const size_t foundNs = numNsInLeftFlank + numNsInRightFlank;

    if (foundNs > maxNsAllowed)
    {
        const string allowed = "Flanks must contain at most " + to_string(maxNsAllowed) + " Ns;";
        const string found = " found " + to_string(foundNs);
        throw runtime_error(allowed + found);
    }

    return leftFlank + structure + rightFlank;
}

static vector<GenomicRegion>
getFeatureLocations(const GraphBlueprint& blueprint, const GraphLocusDecoding& locus, const GenomicRegion& location)
{
    vector<GenomicRegion> locationsOfFlanksAndVariants;
    auto flanks = computeFlanks({ location }, locus.flankLength);
    locationsOfFlanksAndVariants.push_back(flanks.front());
    for (const auto& variant : locus.variants)
    {
        locationsOfFlanksAndVariants.push_back(variant.location);
    }
    locationsOfFlanksAndVariants.push_back(flanks.back());

    vector<GenomicRegion> locationsOfAllFeatures;
    int regionIndex = 0;
    for (const auto& feature : blueprint)
    {
        if (feature.type == GraphBlueprintFeatureType::kInterruption)
        {
            assert(regionIndex != 0 && regionIndex < static_cast<int>(locationsOfFlanksAndVariants.size()));
            const auto& leftRegion = locationsOfFlanksAndVariants[regionIndex];
            const auto& rightRegion = locationsOfFlanksAndVariants[regionIndex + 1];
            locationsOfAllFeatures.emplace_back(leftRegion.contigIndex(), leftRegion.end(), rightRegion.start());
        }
        else
        {
            locationsOfAllFeatures.push_back(locationsOfFlanksAndVariants[regionIndex]);
            ++regionIndex;
        }
    }

    assert(blueprint.size() == locationsOfAllFeatures.size());
    return locationsOfAllFeatures;
}

static NodeLocations
getNodeLocations(const GraphBlueprint& blueprint, const Graph& graph, const vector<GenomicRegion>& featureLocations)
{
    assert(blueprint.size() == featureLocations.size());

    NodeLocations nodeLocations;

    for (int featureIndex = 0; featureIndex != static_cast<int>(blueprint.size()); ++featureIndex)
    {
        const auto& feature = blueprint[featureIndex];
        const auto& featureLocation = featureLocations[featureIndex];

        for (const auto& nodeId : feature.nodeIds)
        {
            const int nodeLength = graph.nodeSeq(nodeId).length();
            GenomicRegion nodeLocation(
                featureLocation.contigIndex(), featureLocation.start(), featureLocation.start() + nodeLength);
            nodeLocations.emplace(std::make_pair(nodeId, nodeLocation));
        }
    }

    return nodeLocations;
}

static GraphBlueprint
getBlueprint(const Reference& reference, const GraphLocusDecoding& locusEncoding, const GenomicRegion& location)
{
    auto locusStructureWithFlanks = addFlanks(reference, locusEncoding.structure, location, locusEncoding.flankLength);
    GraphBlueprint blueprint = decodeFeaturesFromRegex(locusStructureWithFlanks);
    return blueprint;
}

static ReferenceGraph
getGraph(const GraphBlueprint& blueprint, const GraphLocusDecoding& locus, const GenomicRegion& location)
{
    graphtools::Graph graph = makeRegionGraph(blueprint, locus.id);
    auto featureLocations = getFeatureLocations(blueprint, locus, location);
    NodeLocations nodeLocations = getNodeLocations(blueprint, graph, featureLocations);

    return { graph, nodeLocations };
}

static GenotyperParameters getGenotyperParams(const GraphLocusDecoding& encoding)
{
    GenotyperParameters params;
    if (encoding.errorRate)
    {
        params.errorRate = *encoding.errorRate;
    }
    if (encoding.likelihoodRatioThreshold)
    {
        params.likelihoodRatioThreshold = *encoding.likelihoodRatioThreshold;
    }
    if (encoding.minLocusCoverage)
    {
        params.minLocusCoverage = *encoding.minLocusCoverage;
    }

    return params;
}

/*
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
} */

static GraphVariantClassification classifyVariant(
    GraphBlueprintFeatureType featureType, const string& variantTypeFromUser, const GenomicRegion referenceRegion)
{
    if (featureType == GraphBlueprintFeatureType::kInsertionOrDeletion)
    {
        if (referenceRegion.length() == 0)
        {
            return { GraphVariantClassification::Type::kSmallVariant, GraphVariantClassification::Subtype::kInsertion };
        }
        else
        {
            return { GraphVariantClassification::Type::kSmallVariant, GraphVariantClassification::Subtype::kDeletion };
        }
    }
    else if (featureType == GraphBlueprintFeatureType::kSwap)
    {
        if (variantTypeFromUser == "SMN")
        {
            return { GraphVariantClassification::Type::kSmallVariant, GraphVariantClassification::Subtype::kSMN };
        }
        else
        {
            return { GraphVariantClassification::Type::kSmallVariant, GraphVariantClassification::Subtype::kSwap };
        }
    }
    else if (variantTypeFromUser == "Repeat")
    {
        return { GraphVariantClassification::Type::kRepeat, GraphVariantClassification::Subtype::kCommonRepeat };
    }
    else if (variantTypeFromUser == "RareRepeat")
    {
        return { GraphVariantClassification::Type::kRepeat, GraphVariantClassification::Subtype::kRareRepeat };
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

void addVariants(
    const Reference& reference, GraphLocusSpec& locusSpec, const GraphLocusDecoding& locusEncoding,
    const GraphBlueprint& blueprint)
{
    int variantIndex = 0;
    for (const auto& feature : blueprint)
    {
        if (doesFeatureDefineVariant(feature.type))
        {
            auto variantEncoding = locusEncoding.variants.at(variantIndex);
            optional<NodeId> referenceNode = determineReferenceNode(feature, reference, variantEncoding.location);

            // VariantClassification classification(variantType, variantSubtype);

            auto classification = classifyVariant(feature.type, variantEncoding.type, variantEncoding.location);

            locusSpec.addVariant(
                variantEncoding.id, classification, variantEncoding.location, feature.nodeIds, referenceNode);

            ++variantIndex;
        }
    }
}

/*

// TODO: Adapt for new code
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
} */

GraphLocusSpec decode(const Reference& reference, const GraphLocusDecoding& locusEncoding)
{
    auto locusLocation = getLocusLocation(locusEncoding);
    auto copyNumberBySex = getCopyNumber(reference.contigInfo().getContigName(locusLocation.contigIndex()));
    auto analysisRegions = getAnalysisRegions(locusEncoding, locusLocation);
    auto blueprint = getBlueprint(reference, locusEncoding, locusLocation);
    auto graph = getGraph(blueprint, locusEncoding, locusLocation);
    auto genotyperParams = getGenotyperParams(locusEncoding);

    GraphLocusSpec locusSpec(locusEncoding.id, copyNumberBySex, analysisRegions, graph, genotyperParams);
    addVariants(reference, locusSpec, locusEncoding, blueprint);

    return locusSpec;
}
}
