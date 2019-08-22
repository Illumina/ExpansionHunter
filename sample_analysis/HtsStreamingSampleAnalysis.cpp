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

#include "sample_analysis/HtsStreamingSampleAnalysis.hh"

#include <memory>
#include <unordered_map>

#include "common/HtsHelpers.hh"
#include "region_analysis/LocusAnalyzer.hh"
#include "sample_analysis/GenomeQueryCollection.hh"
#include "sample_analysis/HtsFileStreamer.hh"

using graphtools::AlignmentWriter;
using std::map;
using std::string;
using std::unordered_map;
using std::vector;

namespace ehunter
{

SampleFindings htsStreamingSampleAnalysis(
    const InputPaths& inputPaths, Sex sampleSex, const RegionCatalog& regionCatalog, AlignmentWriter& bamletWriter)
{
    vector<std::unique_ptr<LocusAnalyzer>> locusAnalyzers = initializeLocusAnalyzers(regionCatalog, bamletWriter);
    GenomeQueryCollection genomeQuery(locusAnalyzers);

    using ReadCatalog = std::unordered_map<std::string, Read>;
    ReadCatalog unpairedReads;

    htshelpers::HtsFileStreamer readStreamer(inputPaths.htsFile());
    while (readStreamer.trySeekingToNextPrimaryAlignment() && readStreamer.isStreamingAlignedReads())
    {
        const bool isReadNearTargetRegion = genomeQuery.targetRegionMask.query(
            readStreamer.currentReadContigId(), readStreamer.currentReadPosition());
        const bool isMateNearTargetRegion = genomeQuery.targetRegionMask.query(
            readStreamer.currentMateContigId(), readStreamer.currentMatePosition());
        if (!isReadNearTargetRegion && !isMateNearTargetRegion)
        {
            continue;
        }

        LinearAlignmentStats alignmentStats;
        Read read = readStreamer.decodeRead(alignmentStats);
        if (!alignmentStats.isPaired)
        {
            continue;
        }

        const auto mateIterator = unpairedReads.find(read.fragmentId());
        if (mateIterator == unpairedReads.end())
        {
            unpairedReads.emplace(std::make_pair(read.fragmentId(), std::move(read)));
            continue;
        }
        Read mate = std::move(mateIterator->second);
        unpairedReads.erase(mateIterator);

        const int64_t readEnd = readStreamer.currentReadPosition() + read.sequence().length();
        const int64_t mateEnd = readStreamer.currentMatePosition() + mate.sequence().length();

        vector<AnalyzerBundle> analyzerBundles = genomeQuery.analyzerFinder.query(
            readStreamer.currentReadContigId(), readStreamer.currentReadPosition(), readEnd,
            readStreamer.currentMateContigId(), readStreamer.currentMatePosition(), mateEnd);

        for (auto& analyzerBundle : analyzerBundles)
        {
            const auto analyzerPtr = analyzerBundle.locusAnalyzerPtr;

            switch (analyzerBundle.inputType)
            {
            case AnalyzerInputType::kBothReads:
                analyzerPtr->processMates(std::move(read), std::move(mate), analyzerBundle.regionType);
                break;
            case AnalyzerInputType::kReadOnly:
                analyzerPtr->processMates(std::move(read), boost::none, analyzerBundle.regionType);
                break;
            case AnalyzerInputType::kMateOnly:
                analyzerPtr->processMates(std::move(mate), boost::none, analyzerBundle.regionType);
                break;
            }
        }
    }

    SampleFindings sampleFindings;
    for (auto& locusAnalyzer : locusAnalyzers)
    {
        auto locusFindings = locusAnalyzer->analyze(sampleSex);
        sampleFindings.emplace(std::make_pair(locusAnalyzer->locusId(), std::move(locusFindings)));
    }

    return sampleFindings;
}

}
