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

#include "sample_analysis/HtsSeekingSampleAnalyzer.hh"

#include <cassert>
#include <memory>
#include <unordered_map>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>

#include "thirdparty/spdlog/spdlog.h"

#include "reads/ReadPairs.hh"
#include "region_analysis/RegionAnalyzer.hh"
#include "sample_analysis/HtsFileSeeker.hh"
#include "sample_analysis/IndexBasedDepthEstimate.hh"
#include "sample_analysis/MateExtractor.hh"

namespace ehunter
{

using boost::optional;
using graphtools::AlignmentWriter;
using htshelpers::HtsFileSeeker;
using std::ostream;
using std::string;
using std::unordered_map;
using std::vector;

using AlignmentStatsCatalog = unordered_map<ReadId, LinearAlignmentStats, boost::hash<ReadId>>;

static ReadPairs collectReads(
    const vector<GenomicRegion>& regions, AlignmentStatsCatalog& alignmentStatsCatalog, HtsFileSeeker& fileHopper)
{
    ReadPairs readPairs;
    auto console = spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console");

    for (const auto& region : regions)
    {
        int numReadsBeforeCollection = readPairs.NumReads();
        fileHopper.setRegion(region);
        while (fileHopper.trySeekingToNextPrimaryAlignment())
        {
            LinearAlignmentStats alignmentStats;
            Read read = fileHopper.decodeRead(alignmentStats);
            alignmentStatsCatalog.emplace(std::make_pair(read.readId(), alignmentStats));
            readPairs.Add(std::move(read));
        }
        console->debug("Collected {} reads from {}", readPairs.NumReads() - numReadsBeforeCollection, region);
    }
    return readPairs;
}

bool checkIfMatesWereMappedNearby(const LinearAlignmentStats& alignmentStats)
{
    const int maxMateDistance = 1000;
    if ((alignmentStats.chromId == alignmentStats.mateChromId)
        && (std::abs(alignmentStats.pos - alignmentStats.matePos) < maxMateDistance))
    {
        return true;
    }
    return false;
}

void recoverMates(const string& htsFilePath, const AlignmentStatsCatalog& alignmentStatsCatalog, ReadPairs& readPairs)
{
    htshelpers::MateExtractor mateExtractor(htsFilePath);

    for (auto& fragmentIdAndReadPair : readPairs)
    {
        ReadPair& readPair = fragmentIdAndReadPair.second;

        if (readPair.numMatesSet() == 2)
        {
            continue;
        }

        const Read& read = readPair.firstMate != boost::none ? *readPair.firstMate : *readPair.secondMate;

        const auto alignmentStatsIterator = alignmentStatsCatalog.find(read.readId());
        if (alignmentStatsIterator == alignmentStatsCatalog.end())
        {
            throw std::logic_error("Cannot recover mate of uncatalogued read");
        }
        const LinearAlignmentStats& alignmentStats = alignmentStatsIterator->second;

        if (!checkIfMatesWereMappedNearby(alignmentStats))
        {
            optional<Read> optionalMate = mateExtractor.extractMate(read, alignmentStats);
            if (optionalMate != boost::none)
            {
                readPairs.AddMateToExistingRead(*optionalMate);
            }
            else
            {
                auto console = spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console");
                console->warn("Could not recover the mate of {}", read.readId());
            }
        }
    }
}

static RegionFindings analyzeRegion(
    const ReadPairs& readPairs, const ReadPairs& offtargetReadPairs, const LocusSpecification& regionSpec,
    const SampleParameters& sampleParams, const HeuristicParameters& heuristicParams, AlignmentWriter& alignmentWriter)
{
    RegionAnalyzer regionAnalyzer(regionSpec, heuristicParams, alignmentWriter);

    for (const auto fragmentIdAndReads : readPairs)
    {
        const auto& readPair = fragmentIdAndReads.second;
        if (readPair.numMatesSet() == 2)
        {
            regionAnalyzer.processMates(*readPair.firstMate, *readPair.secondMate);
        }
    }

    for (const auto fragmentIdAndReads : offtargetReadPairs)
    {
        const auto& readPair = fragmentIdAndReads.second;
        if (readPair.numMatesSet() == 2)
        {
            regionAnalyzer.processOfftargetMates(*readPair.firstMate, *readPair.secondMate);
        }
    }

    return regionAnalyzer.analyze(sampleParams);
}

SampleFindings htsSeekingSampleAnalysis(
    const InputPaths& inputPaths, SampleParameters& sampleParams, const HeuristicParameters& heuristicParams,
    const RegionCatalog& regionCatalog, AlignmentWriter& alignmentWriter)
{
    auto console = spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console");

    if (!sampleParams.isHaplotypeDepthSet())
    {
        const double depth = estimateDepthFromHtsIndex(inputPaths.htsFile(), sampleParams.readLength());

        const double kMinDepthAllowed = 10.0;
        if (depth < kMinDepthAllowed)
        {
            throw std::invalid_argument("Read depth must be at least " + std::to_string(kMinDepthAllowed));
        }

        sampleParams.setHaplotypeDepth(depth / 2);
        console->info("Depth is set to {}", depth);
    }

    HtsFileSeeker htsFileSeeker(inputPaths.htsFile());

    SampleFindings sampleFindings;
    for (const auto& regionIdAndRegionSpec : regionCatalog)
    {
        const string& regionId = regionIdAndRegionSpec.first;
        const LocusSpecification& regionSpec = regionIdAndRegionSpec.second;

        AlignmentStatsCatalog readAlignmentStats;
        ReadPairs targetReadPairs
            = collectReads(regionSpec.targetReadExtractionRegions(), readAlignmentStats, htsFileSeeker);
        recoverMates(inputPaths.htsFile(), readAlignmentStats, targetReadPairs);
        console->debug("Collected {} read pairs from target regions", targetReadPairs.NumCompletePairs());

        AlignmentStatsCatalog offtargetReadAlignmentStatsCatalog;
        ReadPairs offtargetReadPairs = collectReads(
            regionSpec.offtargetReadExtractionRegions(), offtargetReadAlignmentStatsCatalog, htsFileSeeker);
        recoverMates(inputPaths.htsFile(), offtargetReadAlignmentStatsCatalog, offtargetReadPairs);
        console->debug("Collected {} read pairs from offtarget regions", offtargetReadPairs.NumCompletePairs());

        auto regionFindings = analyzeRegion(
            targetReadPairs, offtargetReadPairs, regionSpec, sampleParams, heuristicParams, alignmentWriter);
        sampleFindings.emplace(std::make_pair(regionId, std::move(regionFindings)));
    }

    return sampleFindings;
}

}
