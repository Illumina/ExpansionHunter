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

#pragma once

#include <memory>
#include <vector>

#include "DepthNormalization.hh"
#include "input/CatalogLoading.hh"
#include "locus_spec/LocusSpec.hh"
#include "output/BamletWriter.hh"
#include "reads/Read.hh"
#include "sample_analysis/ModelFinder.hh"
#include "workflow/LinearModel.hh"
#include "workflow/LocusAnalyzer.hh"
#include "workflow/LocusFindings.hh"
#include "workflow/ReadCountAnalyzer.hh"
#include "workflow/RegionModel.hh"

namespace ehunter
{

class CatalogAnalyzer
{
public:
    CatalogAnalyzer(
        const LocusCatalog& locusCatalog, const std::vector<RegionInfo>& normRegionInfo, BamletWriterPtr bamletWriter);
    void analyze(const MappedRead& read, const MappedRead& mate);
    void analyze(const MappedRead& read);
    void collectResults(
        Sex sampleSex, SampleFindings& sampleFindings, boost::optional<DepthNormalizer> genomeDepthNormalizer);
    DepthNormalizer getGenomeDepthNormalizer();

    const std::vector<std::shared_ptr<RegionModel>>& regionModels() const { return regionModels_; }

private:
    std::vector<std::shared_ptr<LocusAnalyzer>> locusAnalyzers_;
    std::vector<std::shared_ptr<RegionModel>> regionModels_;
    std::unique_ptr<ModelFinder> modelFinder_;
    std::vector<RegionInfo> normRegionInfo_;
    std::shared_ptr<LinearModel> linearModel_;
    std::vector<std::shared_ptr<ReadCountAnalyzer>> readCountAnalyzers_;
};
}
