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

class GraphVariantAlignmentStats
{
public:
    GraphVariantAlignmentStats(double leftBreakpointCoverage, double rightBreakpointCoverage)
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

class GraphVariantAlignmentStatsCalculator
{
public:
    explicit GraphVariantAlignmentStatsCalculator(std::vector<graphtools::NodeId> variantNodes);

    void inspect(const graphtools::GraphAlignment& alignment);
    GraphVariantAlignmentStats getStats(int readLength) const;

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

    std::vector<graphtools::NodeId> variantNodes_;
    graphtools::NodeId firstVariantNode_;
    graphtools::NodeId lastVariantNode_;

    int minMatch_ = 10;
    int numReadsOverlappingLeftBreakpoint_ = 0;
    int numReadsOverlappingRightBreakpoint_ = 0;
};

std::ostream& operator<<(std::ostream& out, const GraphVariantAlignmentStats& stats);

}
