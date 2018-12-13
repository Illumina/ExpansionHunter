//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
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

#include "common/Common.hh"

using std::string;

namespace ehunter
{

Sex decodeSampleSex(const std::string& encoding)
{
    if (encoding == "male")
    {
        return Sex::kMale;
    }
    else if (encoding == "female")
    {
        return Sex::kFemale;
    }
    else
    {
        throw std::invalid_argument(encoding + " is invalid sex; must be either male or female");
    }
}

std::ostream& operator<<(std::ostream& out, ReadType readType)
{
    switch (readType)
    {
    case ReadType::kFlanking:
        out << "FLANKING";
        break;
    case ReadType::kRepeat:
        out << "INREPEAT";
        break;
    case ReadType::kSpanning:
        out << "SPANNING";
        break;
    case ReadType::kOther:
        out << "OTHER";
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, AlleleCount alleleCount)
{
    switch (alleleCount)
    {
    case AlleleCount::kZero:
        out << "Zero alleles";
        break;
    case AlleleCount::kOne:
        out << "One allele";
        break;
    case AlleleCount::kTwo:
        out << "Two alleles";
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, NumericInterval numericInterval)
{
    out << numericInterval.start() << "-" << numericInterval.end();
    return out;
}

}
