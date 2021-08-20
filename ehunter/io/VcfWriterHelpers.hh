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

#include <deque>
#include <set>
#include <string>
#include <utility>

#include "core/Common.hh"
#include "locus/LocusSpecification.hh"
#include "locus/VariantFindings.hh"

namespace ehunter
{

class VcfAlleleFields
{
public:
    VcfAlleleFields(int referenceSize)
        : referenceSize_(referenceSize)
    {
    }

    void addAlleleInfo(
        int alleleSize, ReadType source, NumericInterval confidenceInterval, int spanningReadCount,
        int flankingReadCount, int repeatReadCount);

    std::string encode() const;

private:
    void addRefAlleleInfo(
        int alleleSize, ReadType source, NumericInterval confidenceInterval, int spanningReadCount,
        int flankingReadCount, int repeatReadCount);

    void addAltAlleleInfo(
        int alleleSize, ReadType source, NumericInterval confidenceInterval, int spanningReadCount,
        int flankingReadCount, int repeatReadCount);

    int referenceSize_;
    std::deque<int> genotype_;
    std::deque<ReadType> sources_;
    std::deque<int> alleleSizes_;
    std::deque<NumericInterval> confidenceIntervals_;
    std::deque<int> spanningReadCounts_;
    std::deque<int> flankingReadCounts_;
    std::deque<int> repeatReadCounts_;
};

}
