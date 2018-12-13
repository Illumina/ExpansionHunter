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
using htshelpers::HtsFileSeeker;
using reads::LinearAlignmentStats;
using reads::Read;
using reads::ReadPairs;
using std::ostream;
using std::string;
using std::unordered_map;
using std::vector;

using AlignmentStatsCatalog = unordered_map<string, LinearAlignmentStats>;

static ReadPairs
collectReads(const vector<Region>& regions, AlignmentStatsCatalog& alignmentStatsCatalog, HtsFileSeeker& fileHopper)
{
    ReadPairs readPairs;
    auto console = spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console");
    int preReads = 0;
    for (const auto& region : regions)
    {
        fileHopper.setRegion(region);
        while (fileHopper.trySeekingToNextPrimaryAlignment())
        {
            LinearAlignmentStats alignmentStats;
            Read read = fileHopper.decodeRead(alignmentStats);
            alignmentStatsCatalog.emplace(std::make_pair(read.readId(), alignmentStats));
            readPairs.Add(std::move(read));
        }
        console->info("Collected {} reads from {}", readPairs.NumReads() - preReads, region);
        preReads += readPairs.NumReads();
    }
    return readPairs;
}

bool checkIfMatesWereMappedNearby(const LinearAlignmentStats& alignmentStats)
{
    const int maxMateDistance = 1000;
    if ((alignmentStats.chrom_id == alignmentStats.mate_chrom_id)
        && (std::abs(alignmentStats.pos - alignmentStats.mate_pos) < maxMateDistance))
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
        reads::ReadPair& readPair = fragmentIdAndReadPair.second;

        if (readPair.first_mate.isSet() && readPair.second_mate.isSet())
        {
            continue;
        }

        assert(readPair.first_mate.isSet() || readPair.second_mate.isSet());
        const Read& read = readPair.first_mate.isSet() ? readPair.first_mate : readPair.second_mate;

        const auto alignmentStatsIterator = alignmentStatsCatalog.find(read.readId());
        if (alignmentStatsIterator == alignmentStatsCatalog.end())
        {
            throw std::logic_error("Cannot recover mate of uncatalogued read");
        }
        const LinearAlignmentStats& alignmentStats = alignmentStatsIterator->second;

        if (!checkIfMatesWereMappedNearby(alignmentStats))
        {
            Read mate = mateExtractor.extractMate(read, alignmentStats);
            if (mate.isSet())
            {
                readPairs.AddMateToExistingRead(mate);
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
    const SampleParameters& sampleParams, const HeuristicParameters& heuristicParams, ostream& alignmentStream)
{
    alignmentStream << regionSpec.regionId() << ":" << std::endl;
    RegionAnalyzer regionAnalyzer(regionSpec, sampleParams, heuristicParams, alignmentStream);

    for (const auto fragmentIdAndReads : readPairs)
    {
        const auto& readPair = fragmentIdAndReads.second;
        if (readPair.first_mate.isSet() && readPair.second_mate.isSet())
        {
            regionAnalyzer.processMates(readPair.first_mate, readPair.second_mate);
        }
    }

    for (const auto fragmentIdAndReads : offtargetReadPairs)
    {
        const auto& readPair = fragmentIdAndReads.second;
        if (readPair.first_mate.isSet() && readPair.second_mate.isSet())
        {
            regionAnalyzer.processOfftargetMates(readPair.first_mate, readPair.second_mate);
        }
    }

    return regionAnalyzer.genotype();
}

SampleFindings htsSeekingSampleAnalysis(
    const InputPaths& inputPaths, SampleParameters& sampleParams, const HeuristicParameters& heuristicParams,
    const RegionCatalog& regionCatalog, std::ostream& alignmentStream)
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
        vector<Region> targetRegions;
        const auto& referenceLoci = regionSpec.referenceLoci();
        auto extendRegion = [=](Region region) { return region.extend(heuristicParams.regionExtensionLength()); };
        std::transform(referenceLoci.begin(), referenceLoci.end(), std::back_inserter(targetRegions), extendRegion);
        AlignmentStatsCatalog readAlignmentStats;
        ReadPairs targetReadPairs = collectReads(targetRegions, readAlignmentStats, htsFileSeeker);
        recoverMates(inputPaths.htsFile(), readAlignmentStats, targetReadPairs);
        console->info("Collected {} read pairs from target regions", targetReadPairs.NumCompletePairs());

        AlignmentStatsCatalog offtargetReadAlignmentStatsCatalog;
        ReadPairs offtargetReadPairs
            = collectReads(regionSpec.offtargetLoci(), offtargetReadAlignmentStatsCatalog, htsFileSeeker);
        recoverMates(inputPaths.htsFile(), offtargetReadAlignmentStatsCatalog, offtargetReadPairs);
        console->info("Collected {} read pairs from offtarget regions", offtargetReadPairs.NumCompletePairs());

        auto regionFindings = analyzeRegion(
            targetReadPairs, offtargetReadPairs, regionSpec, sampleParams, heuristicParams, alignmentStream);
        sampleFindings.emplace(std::make_pair(regionId, std::move(regionFindings)));
    }

    return sampleFindings;
}

}
