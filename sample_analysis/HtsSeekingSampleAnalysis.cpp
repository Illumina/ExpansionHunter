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

#include "spdlog/spdlog.h"

#include "DepthNormalization.hh"
#include "common/WorkflowContext.hh"
#include "reads/ReadPairs.hh"
#include "sample_analysis/CatalogAnalyzer.hh"
#include "sample_analysis/HtsFileSeeker.hh"
#include "sample_analysis/MateExtractor.hh"
#include "workflow/LocusAnalyzer.hh"

namespace ehunter
{

using boost::optional;
using graphtools::AlignmentWriter;
using htshelpers::HtsFileSeeker;
using std::dynamic_pointer_cast;
using std::ostream;
using std::shared_ptr;
using std::string;
using std::unique_ptr;
using std::unordered_map;
using std::unordered_set;
using std::vector;

using ReadCatalog = unordered_map<ReadId, MappedRead, boost::hash<ReadId>>;

static bool checkIfMatesWereMappedNearby(const MappedRead& read)
{
    const int kMaxDistance = 1000;
    return (read.contigIndex() == read.mateContigIndex()) && (std::abs(read.pos() - read.matePos()) < kMaxDistance);
}

static void recoverMates(htshelpers::MateExtractor& mateExtractor, ReadPairs& readPairs)
{
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

static int getReadCountCap(const vector<GenomicRegion>& regionsWithReads)
{
    int readCountCap;
    // hardcoded for now
    int sampleDepth = 100;
    int readLength = 150;
    float depthMultiplier = 10;

    int64_t regionLength = 0;
    for (const auto& regionWithReads : regionsWithReads)
    {
        regionLength += regionWithReads.length();
    }

    readCountCap = regionLength / (float)readLength * sampleDepth * depthMultiplier;
    return readCountCap;
}

static ReadPairs collectReads(
    const vector<GenomicRegion>& regionsWithReads, HtsFileSeeker& htsFileSeeker, htshelpers::MateExtractor& mateExtractor)
{
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
    recoverMates(mateExtractor, readPairs);
    const int numReadsAfterRecovery = readPairs.NumReads() - numReadsBeforeRecovery;
    spdlog::debug("Recovered {} reads", numReadsAfterRecovery);

    return readPairs;
}

SampleFindings htsSeekingSampleAnalysis(
    const InputPaths& inputPaths, Sex sampleSex, const LocusCatalog& regionCatalog,
    const vector<RegionInfo>& normRegionInfo, BamletWriterPtr bamletWriter)
{
    HtsFileSeeker htsFileSeeker(inputPaths.htsFile(), inputPaths.reference());
    htshelpers::MateExtractor mateExtractor(inputPaths.htsFile(), inputPaths.reference());
    
	CatalogAnalyzer normRegionAnalyzer({ {} }, normRegionInfo, bamletWriter);
    std::vector<RegionDepthInfo> normDepthInfo;
    for (RegionInfo regionInfo : normRegionInfo)
    {
        ReadPairs readPairs = collectReads({ regionInfo.region }, htsFileSeeker, mateExtractor);

        for (const auto& fragmentIdAndReadPair : readPairs)
        {
            const auto& readPair = fragmentIdAndReadPair.second;
            if (readPair.numMatesSet() == 2)
            {
                const MappedRead& read = *readPair.firstMate;
                const MappedRead& mate = *readPair.secondMate;
                normRegionAnalyzer.analyze(read, mate);
            }
            else
            {
                const MappedRead& read = readPair.firstMate ? *readPair.firstMate : *readPair.secondMate;
                normRegionAnalyzer.analyze(read);
            }
        }
    }

    DepthNormalizer genomeDepthNormalizer = normRegionAnalyzer.getGenomeDepthNormalizer();

    SampleFindings sampleFindings;

    for (const auto& locusIdAndRegionSpec : regionCatalog)
    {
        const auto& locusId = locusIdAndRegionSpec.first;
        const auto& locusSpec = locusIdAndRegionSpec.second;

        ReadPairs readPairs;
        readPairs = collectReads(locusSpec->regionsWithReads(), htsFileSeeker, mateExtractor);

        CatalogAnalyzer catalogAnalyzer({ { locusId, locusSpec } }, std::vector<RegionInfo> {}, bamletWriter);

        for (const auto& fragmentIdAndReadPair : readPairs)
        {
            const auto& readPair = fragmentIdAndReadPair.second;
            if (readPair.numMatesSet() == 2)
            {
                const MappedRead& read = *readPair.firstMate;
                const MappedRead& mate = *readPair.secondMate;
                catalogAnalyzer.analyze(read, mate);
            }
            else
            {
                const MappedRead& read = readPair.firstMate ? *readPair.firstMate : *readPair.secondMate;
                catalogAnalyzer.analyze(read);
            }
        }

        catalogAnalyzer.collectResults(sampleSex, sampleFindings, genomeDepthNormalizer);
    }

    return sampleFindings;
}
}
