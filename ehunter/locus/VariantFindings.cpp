//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "locus/VariantFindings.hh"

using std::string;

namespace ehunter
{

GenotypeFilter operator|(GenotypeFilter left, GenotypeFilter right)
{
    return static_cast<GenotypeFilter>(static_cast<unsigned>(left) | static_cast<unsigned>(right));
}

GenotypeFilter operator&(GenotypeFilter left, GenotypeFilter right)
{
    return static_cast<GenotypeFilter>(static_cast<unsigned>(left) & static_cast<unsigned>(right));
}

std::ostream& operator<<(std::ostream& out, const RepeatFindings& repeatFindings)
{
    out << "Genotype: ";
    if (repeatFindings.optionalGenotype())
    {
        out << *repeatFindings.optionalGenotype();
    }
    else
    {
        out << "N/A";
    }
    out << "; table of spanning reads: " << repeatFindings.countsOfSpanningReads()
        << "; table of flanking reads: " << repeatFindings.countsOfFlankingReads()
        << "; table of inrepeat reads: " << repeatFindings.countsOfInrepeatReads();
    return out;
}

}
