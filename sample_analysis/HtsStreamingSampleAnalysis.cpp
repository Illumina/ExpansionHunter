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
#include <unordered_set>

#include "common/HtsHelpers.hh"
#include "common/WorkflowContext.hh"
#include "sample_analysis/GenomeQueryCollection.hh"
#include "sample_analysis/HtsFileStreamer.hh"
#include "sample_analysis/ReadDispatch.hh"
#include "workflow/LocusAnalyzer.hh"
#include "workflow/WorkflowBuilder.hh"

using graphtools::AlignmentWriter;
using std::map;
using std::shared_ptr;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

namespace ehunter
{

SampleFindings htsStreamingSampleAnalysis(
    const InputPaths& inputPaths, Sex /*sampleSex*/, const RegionCatalog& regionCatalog,
    AlignmentWriter& /*bamletWriter*/)
{
    WorkflowContext context;
    vector<shared_ptr<LocusAnalyzer>> locusAnalyzers;

    for (const auto& locusIdAndLocusSpec : regionCatalog)
    {
        const auto& locusSpec = locusIdAndLocusSpec.second;
        locusAnalyzers.push_back(buildLocusWorkflow(locusSpec, context.heuristics()));
    }

    vector<shared_ptr<RegionModel>> regionModelPtrs = extractRegionModels(locusAnalyzers);

    //= initializeLocusAnalyzers(regionCatalog, bamletWriter);

    GenomeQueryCollection genomeQuery(regionModelPtrs);

    using ReadCatalog = std::unordered_map<std::string, MappedRead>;
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

        MappedRead read = readStreamer.decodeRead();
        if (!read.isPaired())
        {
            continue;
        }

        const auto mateIterator = unpairedReads.find(read.fragmentId());
        if (mateIterator == unpairedReads.end())
        {
            unpairedReads.emplace(std::make_pair(read.fragmentId(), std::move(read)));
            continue;
        }
        MappedRead mate = std::move(mateIterator->second);
        unpairedReads.erase(mateIterator);

        const int64_t readEnd = readStreamer.currentReadPosition() + read.sequence().length();
        const int64_t mateEnd = readStreamer.currentMatePosition() + mate.sequence().length();

        unordered_set<RegionModel*> readModels = genomeQuery.analyzerFinder.query(
            readStreamer.currentReadContigId(), readStreamer.currentReadPosition(), readEnd);

        unordered_set<RegionModel*> mateModels = genomeQuery.analyzerFinder.query(
            readStreamer.currentMateContigId(), readStreamer.currentMateContigId(), mateEnd);

        readModels.insert(mateModels.begin(), mateModels.end());

        for (auto model : readModels)
        {
            model->analyze(read, mate);
        }
    }

    SampleFindings sampleFindings;
    // for (auto& locusAnalyzer : regionModelPtrs)
    //{
    //    auto locusFindings = locusAnalyzer->analyze(sampleSex);
    //    sampleFindings.emplace(std::make_pair(locusAnalyzer->locusId(), std::move(locusFindings)));
    //}

    return sampleFindings;
}
}
