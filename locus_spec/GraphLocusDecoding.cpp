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

#include "locus_spec/GraphBlueprint.hh"
#include "locus_spec/RegionGraph.hh"

using graphtools::Graph;
using std::runtime_error;
using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

static GenomicRegion getLocusLocation(const GraphLocusEncoding& locusEncoding)
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

static AnalysisRegions getAnalysisRegions(const GraphLocusEncoding& encoding, const GenomicRegion& locusLocation)
{
    AnalysisRegions regions;
    regions.regionsWithReads = encoding.regionsWithReads;
    regions.offtargetRegionsWithReads = encoding.offtargetRegionsWithReads;

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
getFeatureLocations(const GraphBlueprint& blueprint, const GraphLocusEncoding& locus, const GenomicRegion& location)
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

static ReferenceGraph
getGraph(const Reference& reference, const GraphLocusEncoding& locus, const GenomicRegion& location)
{
    auto locusStructureWithFlanks = addFlanks(reference, locus.structure, location, locus.flankLength);
    GraphBlueprint blueprint = decodeFeaturesFromRegex(locusStructureWithFlanks);
    graphtools::Graph graph = makeRegionGraph(blueprint, locus.id);
    auto featureLocations = getFeatureLocations(blueprint, locus, location);
    NodeLocations nodeLocations = getNodeLocations(blueprint, graph, featureLocations);

    return { graph, nodeLocations };
}

static GenotyperParameters getGenotyperParams(const GraphLocusEncoding& encoding)
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

GraphLocusSpec decode(const Reference& reference, const GraphLocusEncoding& encoding)
{
    auto locusLocation = getLocusLocation(encoding);
    auto copyNumberBySex = getCopyNumber(reference.contigInfo().getContigName(locusLocation.contigIndex()));
    auto analysisRegions = getAnalysisRegions(encoding, locusLocation);
    auto graph = getGraph(reference, encoding, locusLocation);
    auto genotyperParams = getGenotyperParams(encoding);

    GraphLocusSpec locusSpec(encoding.id, copyNumberBySex, analysisRegions, graph, genotyperParams);
    return locusSpec;
}

}
