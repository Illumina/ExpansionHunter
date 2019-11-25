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
    WorkflowContext context;

    for (const auto& locusIdAndLocusSpec : locusCatalog)
    {
        const auto& locusSpec = locusIdAndLocusSpec.second;
        shared_ptr<CnvLocusSpec> cnvLocusSpec = dynamic_pointer_cast<CnvLocusSpec>(locusSpec);
        shared_ptr<GraphLocusSpec> graphLocusSpec = dynamic_pointer_cast<GraphLocusSpec>(locusSpec);
        shared_ptr<ParalogLocusSpec> paralogLocusSpec = dynamic_pointer_cast<ParalogLocusSpec>(locusSpec);
        if (graphLocusSpec)
        {
            locusAnalyzers_.push_back(buildGraphLocusWorkflow(*graphLocusSpec, context.heuristics(), bamletWriter));
        }
        else if (cnvLocusSpec)
        {
            locusAnalyzers_.push_back(buildCnvLocusWorkflow(*cnvLocusSpec, context.heuristics()));
        }
        else if (paralogLocusSpec)
        {
            locusAnalyzers_.push_back(buildParalogLocusWorkflow(*paralogLocusSpec, context.heuristics()));
        }
    }

    regionModels_ = extractRegionModels(locusAnalyzers_);

    int mapqCutoff = context.heuristics().qualityCutoffForGoodBaseCall();
    int regionExtensionLength = context.heuristics().regionExtensionLength();
    for (RegionInfo regionInfo : normRegionInfo_)
    {
        auto region = regionInfo.region;
        GenomicRegion expandedRegion = GenomicRegion(
            region.contigIndex(), region.start() - regionExtensionLength, region.end() + regionExtensionLength);
        auto linearModel = make_shared<LinearModel>(std::vector<GenomicRegion>{ expandedRegion });
        auto readCounter = make_shared<ReadCounter>(linearModel, std::vector<GenomicRegion>{ region }, mapqCutoff);
        linearModel->addFeature(readCounter.get());
        regionModels_.push_back(linearModel);
        normalizationRegionAnalyzers_.push_back(
            make_shared<ReadCountAnalyzer>(CopyNumberBySex::kTwoInFemaleTwoInMale, readCounter));
    }

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

DepthNormalizer CatalogAnalyzer::getGenomeDepthNormalizer()
{
    std::vector<RegionDepthInfo> normRegionDepthInfo;
    auto regionInfo = normRegionInfo_.begin();
    for (auto normalizationRegionAnalyzer : normalizationRegionAnalyzers_)
    {
        int readCount = normalizationRegionAnalyzer->count();
        double gc = (*regionInfo).gc;
        GenomicRegion region = (*regionInfo).region;
        int regionLength = region.end() - region.start();
        double normDepth = (double)readCount / (double)regionLength;
        normRegionDepthInfo.push_back(RegionDepthInfo(gc, normDepth));
        std::advance(regionInfo, 1);
    }

    return DepthNormalizer(normRegionDepthInfo);
}

void CatalogAnalyzer::collectResults(
    Sex sampleSex, SampleFindings& sampleFindings, boost::optional<DepthNormalizer> genomeDepthNormalizer)
{
    if (!genomeDepthNormalizer)
    {
        genomeDepthNormalizer = getGenomeDepthNormalizer();
    }
    for (auto& locusAnalyzer : locusAnalyzers_)
    {
        auto locusFindings = locusAnalyzer->analyze(sampleSex, genomeDepthNormalizer);
        sampleFindings.emplace(std::make_pair(locusAnalyzer->locusId(), std::move(locusFindings)));
    }
}
}
