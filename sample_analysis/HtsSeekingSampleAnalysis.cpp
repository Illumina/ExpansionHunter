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

#include "sample_analysis/HtsSeekingSampleAnalysis.hh"

#include <cassert>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>

#include "thirdparty/spdlog/spdlog.h"

#include "reads/ReadPairs.hh"
#include "region_analysis/LocusAnalyzer.hh"
#include "sample_analysis/AnalyzerFinder.hh"
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
using std::unique_ptr;
using std::unordered_map;
using std::vector;

namespace
{
    using AlignmentStatsCatalog = unordered_map<ReadId, LinearAlignmentStats, boost::hash<ReadId>>;

    vector<GenomicRegion>
    combineRegions(const vector<GenomicRegion>& targetRegions, const vector<GenomicRegion>& offtargetRegions)
    {
        vector<GenomicRegion> combinedRegions(targetRegions);
        combinedRegions.insert(combinedRegions.end(), offtargetRegions.begin(), offtargetRegions.end());
        return combinedRegions;
    }

    bool checkIfMatesWereMappedNearby(const LinearAlignmentStats& alignmentStats)
    {
        const int kMaxMateDistance = 1000;
        if ((alignmentStats.chromId == alignmentStats.mateChromId)
            && (std::abs(alignmentStats.pos - alignmentStats.matePos) < kMaxMateDistance))
        {
            return true;
        }
        return false;
    }

    void recoverMates(
        const string& htsFilePath, const string& htsReferencePath, AlignmentStatsCatalog& alignmentStatsCatalog,
        ReadPairs& readPairs)
    {
        htshelpers::MateExtractor mateExtractor(htsFilePath, htsReferencePath);

        for (auto& fragmentIdAndReadPair : readPairs)
        {
            ReadPair& readPair = fragmentIdAndReadPair.second;

            if (readPair.numMatesSet() == 2)
            {
                continue;
            }

            const Read& read = readPair.firstMate ? *readPair.firstMate : *readPair.secondMate;

            const auto alignmentStatsIterator = alignmentStatsCatalog.find(read.readId());
            if (alignmentStatsIterator == alignmentStatsCatalog.end())
            {
                throw std::logic_error("Cannot recover mate of uncatalogued read");
            }
            const LinearAlignmentStats& alignmentStats = alignmentStatsIterator->second;

            if (!checkIfMatesWereMappedNearby(alignmentStats))
            {
                LinearAlignmentStats mateStats;
                optional<Read> optionalMate = mateExtractor.extractMate(read, alignmentStats, mateStats);
                if (optionalMate)
                {
                    const Read& mate = *optionalMate;
                    alignmentStatsCatalog.emplace(std::make_pair(mate.readId(), alignmentStats));
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

    ReadPairs collectCandidateReads(
        const vector<GenomicRegion>& targetRegions, const vector<GenomicRegion>& offtargetRegions,
        AlignmentStatsCatalog& alignmentStatsCatalog, const string& htsFilePath, const string& htsReferencePath)
    {
        auto console = spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console");

        vector<GenomicRegion> regionsWithReads = combineRegions(targetRegions, offtargetRegions);
        HtsFileSeeker htsFileSeeker(htsFilePath, htsReferencePath);
        ReadPairs readPairs;

        for (const auto& regionWithReads : regionsWithReads)
        {
            const int numReadsBeforeCollection = readPairs.NumReads();
            htsFileSeeker.setRegion(regionWithReads);
            while (htsFileSeeker.trySeekingToNextPrimaryAlignment())
            {
                LinearAlignmentStats alignmentStats;
                Read read = htsFileSeeker.decodeRead(alignmentStats);
                if (alignmentStats.isPaired)
                {
                    alignmentStatsCatalog.emplace(std::make_pair(read.readId(), alignmentStats));
                    readPairs.Add(std::move(read));
                }
                else
                {
                    console->warn("Skipping {} because it is unpaired", read.readId());
                }
            }
            const int numReadsCollected = readPairs.NumReads() - numReadsBeforeCollection;
            console->debug("Collected {} reads from {}", numReadsCollected, regionWithReads);
        }

        const int numReadsBeforeRecovery = readPairs.NumReads();
        recoverMates(htsFilePath, htsReferencePath, alignmentStatsCatalog, readPairs);
        const int numReadsAfterRecovery = readPairs.NumReads() - numReadsBeforeRecovery;
        console->debug("Recovered {} reads", numReadsAfterRecovery);

        return readPairs;
    }

    void analyzeReadPair(
        AnalyzerFinder& analyzerFinder, const Read& read, const Read& mate, const AlignmentStatsCatalog& alignmentStats)
    {
        const auto readStatsIter = alignmentStats.find(read.readId());
        const auto mateStatsIter = alignmentStats.find(mate.readId());

        if (readStatsIter == alignmentStats.end() || mateStatsIter == alignmentStats.end())
        {
            throw std::logic_error("Could not to find alignment stats for " + read.fragmentId());
        }

        const LinearAlignmentStats& readStats = readStatsIter->second;
        const LinearAlignmentStats& mateStats = mateStatsIter->second;

        const int64_t readEnd = readStats.pos + read.sequence().length();
        const int64_t mateEnd = mateStats.pos + mate.sequence().length();
        vector<AnalyzerBundle> analyzers = analyzerFinder.query(
            readStats.chromId, readStats.pos, readEnd, mateStats.chromId, mateStats.pos, mateEnd);

        if (analyzers.empty())
        {
            return;
        }

        assert(analyzers.size() == 1);
        const AnalyzerBundle& bundle = analyzers.front();

        if (bundle.inputType == AnalyzerInputType::kBothReads)
        {
            bundle.locusAnalyzerPtr->processMates(read, mate, bundle.regionType);
        }
        else if (bundle.inputType == AnalyzerInputType::kReadOnly)
        {
            bundle.locusAnalyzerPtr->processMates(read, boost::none, bundle.regionType);
        }
        else if (bundle.inputType == AnalyzerInputType::kMateOnly)
        {
            bundle.locusAnalyzerPtr->processMates(mate, boost::none, bundle.regionType);
        }
    }

    void analyzeRead(AnalyzerFinder& analyzerFinder, const Read& read, const AlignmentStatsCatalog& alignmentStats)
    {
        const auto readStatsIter = alignmentStats.find(read.readId());

        if (readStatsIter == alignmentStats.end())
        {
            throw std::logic_error("Could not to find alignment stats for " + read.fragmentId());
        }

        const LinearAlignmentStats& readStats = readStatsIter->second;
        const int64_t readEnd = readStats.pos + read.sequence().length();

        vector<AnalyzerBundle> analyzers = analyzerFinder.query(readStats.chromId, readStats.pos, readEnd);

        if (analyzers.empty())
        {
            return;
        }

        assert(analyzers.size() == 1);
        const AnalyzerBundle& bundle = analyzers.front();
        bundle.locusAnalyzerPtr->processMates(read, boost::none, bundle.regionType);
    }

    void processReads(
        const ReadPairs& candidateReadPairs, const AlignmentStatsCatalog& alignmentStats,
        AnalyzerFinder& analyzerFinder)
    {
        for (const auto& fragmentIdAndReads : candidateReadPairs)
        {
            const auto& readPair = fragmentIdAndReads.second;
            if (readPair.numMatesSet() == 2)
            {
                analyzeReadPair(analyzerFinder, *readPair.firstMate, *readPair.secondMate, alignmentStats);
            }
            else
            {
                const Read& read = readPair.firstMate ? *readPair.firstMate : *readPair.secondMate;
                analyzeRead(analyzerFinder, read, alignmentStats);
            }
        }
    }
}

SampleFindings htsSeekingSampleAnalysis(
    const InputPaths& inputPaths, Sex sampleSex, const HeuristicParameters& heuristicParams,
    const RegionCatalog& regionCatalog, AlignmentWriter& alignmentWriter)
{
    auto console = spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console");

    SampleFindings sampleFindings;
    for (const auto& locusIdAndRegionSpec : regionCatalog)
    {
        const string& locusId = locusIdAndRegionSpec.first;
        console->info("Analyzing {}", locusId);
        const LocusSpecification& locusSpec = locusIdAndRegionSpec.second;

        vector<unique_ptr<LocusAnalyzer>> locusAnalyzers;
        locusAnalyzers.emplace_back(new LocusAnalyzer(locusSpec, heuristicParams, alignmentWriter));
        AnalyzerFinder analyzerFinder(locusAnalyzers);

        AlignmentStatsCatalog alignmentStats;
        ReadPairs readPairs = collectCandidateReads(
            locusSpec.targetReadExtractionRegions(), locusSpec.offtargetReadExtractionRegions(), alignmentStats,
            inputPaths.htsFile(), inputPaths.reference());

        processReads(readPairs, alignmentStats, analyzerFinder);

        auto variantFindings = locusAnalyzers.front()->analyze(sampleSex, boost::none);
        sampleFindings.emplace(locusId, std::move(variantFindings));
    }

    return sampleFindings;
}

}
