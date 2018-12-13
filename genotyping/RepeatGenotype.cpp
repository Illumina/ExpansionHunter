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
