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

#include <boost/optional.hpp>
#include <cstdint>

#include "workflow/LinearFeature.hh"
#include "locus_spec/ParalogLocusSpec.hh"

namespace ehunter
{

class LinearSmallVariant : public LinearFeature
{
public:
    LinearSmallVariant(
        std::shared_ptr<LinearModel> model, SmallVariantLocations locations, SmallVariantBases bases, boost::optional<int> mapqCutoff);
    ~LinearSmallVariant() override = default;
    std::shared_ptr<RegionModel> model() override;
    void summarize(const MappedRead& read, const MappedRead& mate) override;
    void summarize(const MappedRead& read) override;

    std::int64_t numGeneAReads() const { return numGeneAReads_; }
    std::int64_t numGeneBReads() const { return numGeneBReads_; }
    

private:
    std::shared_ptr<LinearModel> model_;
    SmallVariantLocations locations_;
    SmallVariantBases bases_;
    boost::optional<int> mapqCutoff_;

    std::int64_t numGeneAReads_ = 0;
    std::int64_t numGeneBReads_ = 0;
};
}
