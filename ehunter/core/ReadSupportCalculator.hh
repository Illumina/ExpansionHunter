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

#include <cassert>

#include "core/CountTable.hh"

namespace ehunter
{

// Determines counts of informative reads consistent with a given repeat length
class ReadSupportCalculator
{
public:
    ReadSupportCalculator(
        const CountTable& spanningReadCounts, const CountTable& flankingReadCounts,
        const CountTable& inrepeatReadCounts)
        : spanningReadCounts_(spanningReadCounts)
        , flankingReadCounts_(flankingReadCounts)
        , inrepeatReadCounts_(inrepeatReadCounts)

    {
    }

    // A spanning read is consistent with the given repeat allele if it spans same number of repeat units
    int getCountOfConsistentSpanningReads(int haplotypeSize) const;
    // A flanking read is consistent with the given repeat allele if it spans the same or fewer number of repeat units
    int getCountOfConsistentFlankingReads(int haplotypeSize) const;
    // Reports the number of in-repeat reads if the repeat allele is longer than the read length
    int getCountOfConsistentRepeatReads(int haplotypeSize) const;

private:
    const CountTable& spanningReadCounts_;
    const CountTable& flankingReadCounts_;
    const CountTable& inrepeatReadCounts_;
};

}
