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

#pragma once

#include <iostream>
#include <list>
#include <set>

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

namespace ehunter
{

class StrAlignmentStats
{
public:
    StrAlignmentStats(int numReadsOverlappingLeftBreakpoint, int numReadsOverlappingRightBreakpoint)
        : numReadsOverlappingLeftBreakpoint_(numReadsOverlappingLeftBreakpoint)
        , numReadsOverlappingRightBreakpoint_(numReadsOverlappingRightBreakpoint)
    {
    }

    int numReadsOverlappingLeftBreakpoint() const { return numReadsOverlappingLeftBreakpoint_; }
    int numReadsOverlappingRightBreakpoint() const { return numReadsOverlappingRightBreakpoint_; }

private:
    int numReadsOverlappingLeftBreakpoint_;
    int numReadsOverlappingRightBreakpoint_;
};

class StrAlignmentStatsCalculator
{
public:
    explicit StrAlignmentStatsCalculator(graphtools::NodeId strNode);

    void inspect(const std::list<graphtools::GraphAlignment>& alignments);
    StrAlignmentStats getStats() const;

private:
    enum class Flank
    {
        kLeft,
        kRight,
        kBoth,
        kNeither
    };

    Flank classify(const graphtools::GraphAlignment& alignment);

    graphtools::NodeId strNode_;
    int minMatch_ = 10;
    int numReadsOverlappingLeftBreakpoint_ = 0;
    int numReadsOverlappingRightBreakpoint_ = 0;
};

std::ostream& operator<<(std::ostream& out, const StrAlignmentStats& stats);

}
