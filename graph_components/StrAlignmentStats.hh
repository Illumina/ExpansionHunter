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
    StrAlignmentStats(double leftBreakpointCoverage, double rightBreakpointCoverage)
        : leftBreakpointCoverage_(leftBreakpointCoverage)
        , rightBreakpointCoverage_(rightBreakpointCoverage)
    {
    }

    double leftBreakpointCoverage() const { return leftBreakpointCoverage_; }
    double rightBreakpointCoverage() const { return rightBreakpointCoverage_; }

private:
    double leftBreakpointCoverage_;
    double rightBreakpointCoverage_;
};

class StrAlignmentStatsCalculator
{
public:
    explicit StrAlignmentStatsCalculator(graphtools::NodeId strNode);

    void inspect(const std::list<graphtools::GraphAlignment>& alignments);
    StrAlignmentStats getStats(int readLength) const;

private:
    enum class Flank
    {
        kLeft,
        kRight,
        kBoth,
        kNeither
    };

    Flank classify(const graphtools::GraphAlignment& alignment) const;
    double computeBreakpointCoverage(int numReads, int readLength) const;

    graphtools::NodeId strNode_;
    int minMatch_ = 10;
    int numReadsOverlappingLeftBreakpoint_ = 0;
    int numReadsOverlappingRightBreakpoint_ = 0;
};

std::ostream& operator<<(std::ostream& out, const StrAlignmentStats& stats);

}
