//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#pragma once

#include <deque>
#include <set>
#include <string>
#include <utility>

#include "common/common.h"
#include "region_analysis/RepeatFindings.hh"
#include "region_spec/RegionSpec.hh"

class VcfSampleFields
{
public:
    VcfSampleFields(int referenceSize)
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

std::set<int> computeAltRepeatSizes(const RegionCatalog& regionSpecs, const SampleFindings& sampleFindings);
