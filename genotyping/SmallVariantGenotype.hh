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

#include <iostream>
#include <vector>

namespace ehunter
{

enum class AlleleType
{
    kRef,
    kAlt
};

class SmallVariantGenotype
{
public:
    SmallVariantGenotype(AlleleType firstAlleleType, AlleleType secondAlleleType)
    {
        alleleTypes_ = { firstAlleleType, secondAlleleType };
    }
    SmallVariantGenotype(AlleleType alleleType) { alleleTypes_ = { alleleType }; }

    int numAlleles() const { return alleleTypes_.size(); }
    AlleleType firstAlleleType() const { return alleleTypes_.front(); }
    AlleleType secondAlleleType() const { return alleleTypes_.back(); }

    bool isHomRef() const { return firstAlleleType() == AlleleType::kRef && secondAlleleType() == AlleleType::kRef; }
    bool isHomAlt() const { return firstAlleleType() == AlleleType::kAlt && secondAlleleType() == AlleleType::kAlt; }

    bool operator==(const SmallVariantGenotype& other) const { return alleleTypes_ == other.alleleTypes_; }

private:
    std::vector<AlleleType> alleleTypes_;
};

std::ostream& operator<<(std::ostream& out, const SmallVariantGenotype& genotype);

}
