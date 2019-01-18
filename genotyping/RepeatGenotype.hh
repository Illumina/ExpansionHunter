//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
// Concept: Michael Eberle <meberle@illumina.com>
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

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "common/Common.hh"

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
