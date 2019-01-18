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

#include "genotyping/SmallVariantGenotype.hh"

using std::ostream;

namespace ehunter
{

ostream& operator<<(ostream& out, const SmallVariantGenotype& genotype)
{

    if (genotype.numAlleles() == 1)
    {
        out << (genotype.isHomRef() ? "0" : "1");
    }
    else
    {
        if (genotype.isHomRef())
        {
            out << "0/0";
        }
        else if (genotype.isHomAlt())
        {
            out << "1/1";
        }
        else
        {
            out << "0/1";
        }
    }

    return out;
}
}
