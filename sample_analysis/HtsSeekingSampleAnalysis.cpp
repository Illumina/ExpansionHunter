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
#include <unordered_set>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>

#include "thirdparty/spdlog/spdlog.h"

#include "common/WorkflowContext.hh"
#include "reads/ReadPairs.hh"
#include "region/LocusAnalyzer.hh"
#include "region/WorkflowBuilder.hh"
#include "sample_analysis/HtsFileSeeker.hh"
#include "sample_analysis/IndexBasedDepthEstimate.hh"
#include "sample_analysis/MateExtractor.hh"
#include "sample_analysis/ModelFinder.hh"
#include "sample_analysis/ReadDispatch.hh"

namespace ehunter
{

using boost::optional;
using graphtools::AlignmentWriter;
using htshelpers::HtsFileSeeker;
using std::ostream;
using std::shared_ptr;
using std::string;
using std::unique_ptr;
using std::unordered_map;
using std::unordered_set;
using std::vector;

namespace
{
    using ReadCatalog = unordered_map<ReadId, MappedRead, boost::hash<ReadId>>;

    vector<GenomicRegion>
    combineRegions(const vector<GenomicRegion>& targetRegions, const vector<GenomicRegion>& offtargetRegions)
    {
        vector<GenomicRegion> combinedRegions(targetRegions);
        combinedRegions.insert(combinedRegions.end(), offtargetRegions.begin(), offtargetRegions.end());
        return combinedRegions;
    }

    bool checkIfMatesWereMappedNearby(const MappedRead& read)
    {
        const int kMaxDistance = 1000;
        return (read.contigIndex() == read.mateContigIndex()) && (std::abs(read.pos() - read.matePos()) < kMaxDistance);
    }

    void recoverMates(const string& htsFilePath, const string& htsReferencePath, ReadPairs& readPairs)
    {
        htshelpers::MateExtractor mateExtractor(htsFilePath, htsReferencePath);

        for (auto& fragmentIdAndReadPair : readPairs)
        {
            ReadPair& readPair = fragmentIdAndReadPair.second;

            if (readPair.numMatesSet() == 2)
            {
                continue;
            }

            const MappedRead& read = readPair.firstMate ? *readPair.firstMate : *readPair.secondMate;

            if (!checkIfMatesWereMappedNearby(read))
            {
                optional<MappedRead> optionalMate = mateExtractor.extractMate(read);
                if (optionalMate)
                {
                    const MappedRead& mate = *optionalMate;
                    readPairs.AddMateToExistingRead(mate);
                }
                else
                {
                    // TODO: Uncomment
                    // auto console = spdlog::get("console") ? spdlog::get("console") :
                    // spdlog::stderr_color_mt("console"); console->warn("Could not recover the mate of {}",
                    // read.readId());
                }
            }
        }
    }

    int32_t getReadCountCap(vector<GenomicRegion>& regionsWithReads)
    {
        int32_t readCountCap;
        // hardcoded for now
        int32_t sampleDepth = 100;
        int32_t readLength = 150;
        float depthMultiplier = 10;

        int64_t regionLength = 0;
        for (const auto& regionWithReads : regionsWithReads)
        {
            regionLength += regionWithReads.length();
        }

        readCountCap = regionLength / (float)readLength * sampleDepth * depthMultiplier;
        return readCountCap;
    }

    ReadPairs collectCandidateReads(
        const vector<GenomicRegion>& targetRegions, const vector<GenomicRegion>& offtargetRegions,
        const string& htsFilePath, const string& htsReferencePath)
    {
        auto console = spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console");

        vector<GenomicRegion> regionsWithReads = combineRegions(targetRegions, offtargetRegions);
        HtsFileSeeker htsFileSeeker(htsFilePath, htsReferencePath);
        ReadPairs readPairs;

        for (const auto& regionWithReads : regionsWithReads)
        {
            // const int numReadsBeforeCollection = readPairs.NumReads();
            htsFileSeeker.setRegion(regionWithReads);
            while (htsFileSeeker.trySeekingToNextPrimaryAlignment())
            {
                MappedRead read = htsFileSeeker.decodeRead();
                if (read.isPaired())
                {
                    readPairs.Add(std::move(read));
                }
                else
                {
                    // TODO: Uncomment
                    // console->warn("Skipping {} because it is unpaired", read.readId());
                }
            }
            // const int numReadsCollected = readPairs.NumReads() - numReadsBeforeCollection;
            // console->debug("Collected {} reads from {}", numReadsCollected, regionWithReads);
        }

        // add a cap for reads
        if (readPairs.NumReads() > getReadCountCap(regionsWithReads))
        {
            readPairs.Clear();
        }

        const int numReadsBeforeRecovery = readPairs.NumReads();
        recoverMates(htsFilePath, htsReferencePath, readPairs);
        const int numReadsAfterRecovery = readPairs.NumReads() - numReadsBeforeRecovery;
        console->debug("Recovered {} reads", numReadsAfterRecovery);

        return readPairs;
    }

    void analyzeReadPair(ModelFinder& analyzerFinder, const MappedRead& read, const MappedRead& mate)
    {

        const int64_t readEnd = read.pos() + read.sequence().length();
        const int64_t mateEnd = mate.pos() + mate.sequence().length();

        unordered_set<RegionModel*> readModels = analyzerFinder.query(read.contigIndex(), read.pos(), readEnd);
        unordered_set<RegionModel*> mateModels = analyzerFinder.query(mate.contigIndex(), mate.pos(), mateEnd);

        readModels.insert(mateModels.begin(), mateModels.end());
        dispatch(read, mate, readModels);

        /*
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
            bundle.regionPtr->analyze(read, mate);
        }
        else if (bundle.inputType == AnalyzerInputType::kReadOnly)
        {
            bundle.regionPtr->analyze(read, boost::none);
        }
        else if (bundle.inputType == AnalyzerInputType::kMateOnly)
        {
            bundle.regionPtr->analyze(mate, boost::none);
        } */
    }

    /*
    void analyzeRead(ModelFinder& analyzerFinder, const Read& read, const AlignmentStatsCatalog& alignmentStats)
    {
        const auto readStatsIter = alignmentStats.find(read.readId());

        if (readStatsIter == alignmentStats.end())
        {
            throw std::logic_error("Could not to find alignment stats for " + read.fragmentId());
        }

        const LinearAlignmentStats& readStats = readStatsIter->second;
        const int64_t readEnd = readStats.pos + read.sequence().length();

        vector<RegionModel*> analyzers = analyzerFinder.query(readStats.chromId, readStats.pos, readEnd);

        if (analyzers.empty())
        {
            return;
        }

        assert(analyzers.size() == 1);
        const AnalyzerBundle& bundle = analyzers.front();
        bundle.regionPtr->analyze(read, boost::none);
    }
    */

    void processReads(const ReadPairs& candidateReadPairs, ModelFinder& analyzerFinder)
    {
        for (const auto& fragmentIdAndReads : candidateReadPairs)
        {
            const auto& readPair = fragmentIdAndReads.second;
            if (readPair.numMatesSet() == 2)
            {
                analyzeReadPair(analyzerFinder, *readPair.firstMate, *readPair.secondMate);
            }
            // else
            //{
            //    const Read& read = readPair.firstMate ? *readPair.firstMate : *readPair.secondMate;
            //    analyzeRead(analyzerFinder, read, alignmentStats);
            //}
        }
    }
}

SampleFindings htsSeekingSampleAnalysis(
    const InputPaths& inputPaths, Sex /*sampleSex*/, const RegionCatalog& regionCatalog,
    AlignmentWriter& /*alignmentWriter*/)
{
    auto console = spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console");

    SampleFindings sampleFindings;
    for (const auto& locusIdAndRegionSpec : regionCatalog)
    {
        // const string& locusId = locusIdAndRegionSpec.first;
        const LocusSpecification& locusSpec = locusIdAndRegionSpec.second;

        // TODO: Initialize regions
        WorkflowContext context;
        shared_ptr<LocusAnalyzer> locusAnalyzerPtr = buildLocusWorkflow(locusSpec, context.heuristics());

        // TODO: For each region: collect reads and pass them into regions

        // vector<unique_ptr<LocusAnalyzer>> locusAnalyzers;
        // locusAnalyzers.emplace_back(new LocusAnalyzer(locusSpec, alignmentWriter));
        vector<shared_ptr<RegionModel>> regionModelPtrs = extractRegionModels({ locusAnalyzerPtr });
        ModelFinder analyzerFinder(regionModelPtrs);

        ReadPairs readPairs = collectCandidateReads(
            locusSpec.targetReadExtractionRegions(), locusSpec.offtargetReadExtractionRegions(), inputPaths.htsFile(),
            inputPaths.reference());

        processReads(readPairs, analyzerFinder);

        // auto variantFindings = locusAnalyzers.front()->analyze(sampleSex);
        // sampleFindings.emplace(locusId, std::move(variantFindings));
    }

    return sampleFindings;
}
}
