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

#include <cassert>

#include "common/count_table.h"

// Determines counts of informative reads consistent with a given repeat length
class ReadSupportCalculator
{
public:
    ReadSupportCalculator(
        int maxUnitsInRead, const CountTable& spanningReadCounts, const CountTable& flankingReadCounts)
        : maxUnitsInRead_(maxUnitsInRead)
        , spanningReadCounts_(spanningReadCounts)
        , flankingReadCounts_(flankingReadCounts)

    {
        assert(maxUnitsInRead_ > 0);
    }

    // A spanning read is consistent with the given repeat allele if it spans same number of repeat units
    int getCountOfConsistentSpanningReads(int haplotypeSize) const;
    // A flanking read is consistent with the given repeat allele if it spans the same or fewer number of repeat units
    int getCountOfConsistentFlankingReads(int haplotypeSize) const;
    // Reports the number of in-repeat reads if the repeat allele is longer than the read length
    int getCountOfConsistentRepeatReads(int haplotypeSize) const;

private:
    const int maxUnitsInRead_;
    const CountTable& spanningReadCounts_;
    const CountTable& flankingReadCounts_;
};