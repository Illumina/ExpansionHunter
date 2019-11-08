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

#include "sample_analysis/CatalogAnalyzer.hh"
#include "common/WorkflowContext.hh"
#include "workflow/RegionModel.hh"
#include "workflow/WorkflowBuilder.hh"

using std::dynamic_pointer_cast;
using std::shared_ptr;
using std::unordered_set;
using std::vector;

namespace ehunter
{

CatalogAnalyzer::CatalogAnalyzer(
    const LocusCatalog& locusCatalog, const std::vector<RegionInfo>& normRegionInfo, BamletWriterPtr bamletWriter)
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

    WorkflowContext context;

    for (const auto& locusIdAndLocusSpec : locusCatalog)
    {
        const auto& locusSpec = locusIdAndLocusSpec.second;
        shared_ptr<CnvLocusSpec> cnvLocusSpec = dynamic_pointer_cast<CnvLocusSpec>(locusSpec);
        shared_ptr<GraphLocusSpec> graphLocusSpec = dynamic_pointer_cast<GraphLocusSpec>(locusSpec);
        if (graphLocusSpec)
        {
            locusAnalyzers_.push_back(buildGraphLocusWorkflow(*graphLocusSpec, context.heuristics(), bamletWriter));
        }
        else if (cnvLocusSpec)
        {
            locusAnalyzers_.push_back(buildCnvLocusWorkflow(*cnvLocusSpec, genomeDepthNormalizer));
        }
    }

    regionModels_ = extractRegionModels(locusAnalyzers_);
    auto kk = std::unique_ptr<RegionModel>(linearModel_);
    regionModels.emplace_back();

    modelFinder_.reset(new ModelFinder(regionModels_));
}

void CatalogAnalyzer::analyze(const MappedRead& read, const MappedRead& mate)
{
    unordered_set<RegionModel*> readModels = modelFinder_->query(read.contigIndex(), read.pos(), read.approximateEnd());
    unordered_set<RegionModel*> mateModels = modelFinder_->query(mate.contigIndex(), mate.pos(), mate.approximateEnd());

    readModels.insert(mateModels.begin(), mateModels.end());

    for (auto model : readModels)
    {
        model->analyze(read, mate);
    }
}

void CatalogAnalyzer::analyze(const MappedRead& read)
{
    unordered_set<RegionModel*> readModels = modelFinder_->query(read.contigIndex(), read.pos(), read.approximateEnd());
    for (auto model : readModels)
    {
        model->analyze(read);
    }
}

std::vector<RegionDepthInfo> CatalogAnalyzer::summarize()
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

void CatalogAnalyzer::collectResults(Sex sampleSex, SampleFindings& sampleFindings)
{
    for (auto& locusAnalyzer : locusAnalyzers_)
    {
        auto locusFindings = locusAnalyzer->analyze(sampleSex);
        sampleFindings.emplace(std::make_pair(locusAnalyzer->locusId(), std::move(locusFindings)));
    }
}
}
