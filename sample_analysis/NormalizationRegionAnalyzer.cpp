//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Xiao Chen <xchen2@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "NormalizationRegionAnalyzer.hh"

using boost::optional;
using std::unordered_set;
using std::vector;
#include <string>
#include <vector>

namespace ehunter
{

NormalizationRegionAnalyzer::NormalizationRegionAnalyzer(const std::vector<RegionInfo>& normRegionInfo)
    : normRegionInfo_(normRegionInfo)
{
    std::vector<GenomicRegion> normRegions;
    for (RegionInfo regionInfo : normRegionInfo_)
    {
        normRegions.push_back(regionInfo.region);
    }
    linearModel_ = make_shared<LinearModel>(normRegions);
    for (RegionInfo regionInfo : normRegionInfo_)
    {
        std::vector<GenomicRegion> countingRegion = std::vector<GenomicRegion> { regionInfo.region };
        auto readCounter = make_shared<ReadCounter>(linearModel_, countingRegion);
        linearModel_->addFeature(readCounter.get());
        readCountAnalyzers_.push_back(
            make_shared<ReadCountAnalyzer>(CopyNumberBySex::kTwoInFemaleTwoInMale, readCounter));
    }
}

void NormalizationRegionAnalyzer::analyze(const MappedRead& read, const MappedRead& mate)
{
    linearModel_->analyze(read, mate);
}

void NormalizationRegionAnalyzer::analyze(const MappedRead& read) { linearModel_->analyze(read); }

std::vector<RegionDepthInfo> NormalizationRegionAnalyzer::summarize()
{
    std::vector<RegionDepthInfo> normRegionDepthInfo;
    auto regionInfo = normRegionInfo_.begin();
    for (auto readCountAnalyzer : readCountAnalyzers_)
    {
        int readCount = readCountAnalyzer->count();
        double gc = (*regionInfo).gc;
        GenomicRegion region = (*regionInfo).region;
        int regionLength = region.end() - region.start();
        double normDepth = (double)readCount / (double)regionLength;
        normRegionDepthInfo.push_back(RegionDepthInfo(gc, normDepth));
        std::advance(regionInfo, 1);
    }

    return normRegionDepthInfo;
}
}
