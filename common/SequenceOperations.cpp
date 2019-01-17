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

#include "common/SequenceOperations.hh"

using std::string;

namespace ehunter
{

string lowercaseLowQualityBases(const string& bases, const string& quals, int low_base_quality_cutoff)
{
    string cased_bases = bases;
    for (size_t index = 0; index != bases.size(); ++index)
    {
        if (quals[index] - 33 <= low_base_quality_cutoff)
        {
            cased_bases[index] = std::tolower(bases[index]);
        }
    }
    return cased_bases;
}

}
