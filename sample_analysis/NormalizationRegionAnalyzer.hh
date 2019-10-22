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
#pragma once
#include "DepthNormalization.hh"
#include "common/Common.hh"
#include "common/GenomicRegion.hh"
#include "input/CatalogLoading.hh"
#include "reads/Read.hh"
#include "workflow/LinearModel.hh"
#include "workflow/ReadCountAnalyzer.hh"
#include <boost/optional.hpp>
#include <memory>
#include <string>
#include <vector>

namespace ehunter
{

class NormalizationRegionAnalyzer
{
public:
    NormalizationRegionAnalyzer(const std::vector<RegionInfo>& normRegionInfo);
    void analyze(const MappedRead& read, const MappedRead& mate);
    void analyze(const MappedRead& read);
    std::vector<RegionDepthInfo> summarize();

private:
    std::shared_ptr<LinearModel> linearModel_;
    std::vector<std::shared_ptr<ReadCountAnalyzer>> readCountAnalyzers_;
    std::vector<RegionInfo> normRegionInfo_;
};
}