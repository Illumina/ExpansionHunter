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

#include <cstdint>

#include "workflow/LinearFeature.hh"

namespace ehunter
{

class ReadCounter : public LinearFeature
{
public:
    ReadCounter(std::shared_ptr<LinearModel> model, std::vector<GenomicRegion> targetRegions);
    ~ReadCounter() override = default;
    std::shared_ptr<RegionModel> model() override;
    void summarize(MappedRead read, MappedRead mate) override;
    void summarize(MappedRead read) override;

    std::int64_t numReads() const { return numReads_; }
    int getReadLength() const;
    double getDepth() const;

private:
    std::shared_ptr<LinearModel> model_;
    std::vector<GenomicRegion> targetRegions_;

    std::int64_t numReads_ = 0;
    std::int64_t totalReadLength_ = 0;
};

}
