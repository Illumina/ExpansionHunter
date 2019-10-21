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
#include "DepthNormalization.hh"
#include "input/CatalogLoading.hh"
#include "workflow/LinearModel.hh"
#include "workflow/ReadCountAnalyzer.hh"
#include "workflow/RegionModel.hh"

using boost::optional;
using std::unordered_set;
using std::vector;
#include <string>
#include <vector>

namespace ehunter
{

NormalizationRegionAnalyzer::NormalizationRegionAnalyzer(GenomicRegion region)
: region_(region)
{
    std::vector<GenomicRegion> countingRegion = std::vector<GenomicRegion>{region_};
    auto linearModel = make_shared<LinearModel>(countingRegion);
    readCounter_ = make_shared<ReadCounter>(linearModel, countingRegion);
    linearModel->addFeature(readCounter_.get());
    readCountAnalyzer_ = make_shared<ReadCountAnalyzer>(ContigCopyNumber::kTwoInFemaleTwoInMale, readCounter_);
}

void NormalizationRegionAnalyzer::analyze(const MappedRead& read, const MappedRead& mate)
{
    auto model = readCounter_->model();
    model->analyze(read, mate);
    
}

void NormalizationRegionAnalyzer::analyze(const MappedRead& read)
{
    auto model = readCounter_->model();
    model->analyze(read);
}

double NormalizationRegionAnalyzer::summarize()
{
    int readCount = readCountAnalyzer_->count();
    int regionLength = region_.end() - region_.start();
    return (double)readCount / (double)regionLength;
}
}

