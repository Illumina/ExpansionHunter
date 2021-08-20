//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
// Concept: Michael Eberle <meberle@illumina.com>
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

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/Common.hh"

namespace ehunter
{

class RepeatGenotype
{
public:
    RepeatGenotype(int32_t repeatUnitLen, const std::vector<int32_t>& alleleSizes)
        : repeatUnitLen_(repeatUnitLen)
    {
        for (int32_t size : alleleSizes)
        {
            RepeatAllele repeatAllele;
            repeatAllele.numRepeatUnits = size;
            repeatAllele.ci = NumericInterval(size, size);
            alleles_.push_back(repeatAllele);
        }

        assertValidity();
    }

    int32_t repeatUnitLen() const { return repeatUnitLen_; }
    int32_t numAlleles() const { return alleles_.size(); }
    bool isHomozygous() const { return shortAlleleSizeInUnits() == longAlleleSizeInUnits(); }
    int32_t shortAlleleSizeInUnits() const { return alleles_.front().numRepeatUnits; }
    void setShortAlleleSizeInUnits(int32_t numUnits)
    {
        alleles_.front().numRepeatUnits = numUnits;
        alleles_.front().ci = NumericInterval(numUnits, numUnits);
        assertValidity();
    }
    int32_t shortAlleleSizeInBp() const { return alleles_.front().numRepeatUnits * repeatUnitLen_; }
    NumericInterval shortAlleleSizeInUnitsCi() const { return alleles_.front().ci; }
    void setShortAlleleSizeInUnitsCi(int32_t lowerBound, int32_t upperBound)
    {
        const int shortAlleleSize = alleles_.front().numRepeatUnits;
        lowerBound = std::min(lowerBound, shortAlleleSize);
        upperBound = std::max(upperBound, shortAlleleSize);
        alleles_.front().ci = NumericInterval(lowerBound, upperBound);
        assertValidity();
    }

    int32_t longAlleleSizeInUnits() const { return alleles_.back().numRepeatUnits; }
    void setLongAlleleSizeInUnits(int32_t numUnits)
    {
        alleles_.back().numRepeatUnits = numUnits;
        alleles_.back().ci = NumericInterval(numUnits, numUnits);
        assertValidity();
    }

    int32_t longAlleleSizeInBp() const { return alleles_.back().numRepeatUnits * repeatUnitLen_; }
    NumericInterval longAlleleSizeInUnitsCi() const { return alleles_.back().ci; }
    void setLongAlleleSizeInUnitsCi(int32_t lowerBound, int32_t upperBound)
    {
        const int longAlleleSize = alleles_.back().numRepeatUnits;
        lowerBound = std::min(lowerBound, longAlleleSize);
        upperBound = std::max(upperBound, longAlleleSize);
        alleles_.back().ci = NumericInterval(lowerBound, upperBound);
        assertValidity();
    }

    bool operator==(const RepeatGenotype& other) const
    {
        return repeatUnitLen_ == other.repeatUnitLen_ && alleles_ == other.alleles_;
    }

private:
    struct RepeatAllele
    {
        int32_t numRepeatUnits = 0;
        NumericInterval ci;
        bool operator==(const RepeatAllele& other) const
        {
            return numRepeatUnits == other.numRepeatUnits && ci == other.ci;
        }
    };
    void assertValidity() const;
    int32_t repeatUnitLen_;
    std::vector<RepeatAllele> alleles_;
};

std::ostream& operator<<(std::ostream& out, const RepeatGenotype& genotype);

}
