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
