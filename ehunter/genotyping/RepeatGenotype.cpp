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

#include "genotyping/RepeatGenotype.hh"

#include <stdexcept>

using std::ostream;
using std::string;
using std::to_string;

namespace ehunter
{

void RepeatGenotype::assertValidity() const
{
    if (alleles_.empty() || alleles_.size() > 2)
    {
        throw std::logic_error(std::to_string(alleles_.size()) + " is not a valid number of alleles");
    }

    if (shortAlleleSizeInBp() > longAlleleSizeInBp())
    {
        throw std::logic_error("Allele sizes are not ordered");
    }

    for (const RepeatAllele& allele : alleles_)
    {
        const bool isCiOrdered = allele.ci.start() <= allele.ci.end();
        const bool isRepeatSizeInsizeCi
            = allele.ci.start() <= allele.numRepeatUnits && allele.numRepeatUnits <= allele.ci.end();

        if (!isCiOrdered || !isRepeatSizeInsizeCi)
        {
            string ciEncoding = "(" + to_string(allele.ci.start()) + ", " + to_string(allele.ci.end()) + ")";
            string repeatSizeEncoding = to_string(allele.numRepeatUnits);
            throw std::logic_error(ciEncoding + " is invalid CI for repeat of size " + repeatSizeEncoding);
        }
    }
}

ostream& operator<<(ostream& out, const RepeatGenotype& genotype)
{
    const auto& shortAlleleCi = genotype.shortAlleleSizeInUnitsCi();

    out << shortAlleleCi;

    if (genotype.numAlleles() == 2)
    {
        const auto& longAlleleCi = genotype.longAlleleSizeInUnitsCi();

        out << "/" << longAlleleCi;
    }

    return out;
}

}
