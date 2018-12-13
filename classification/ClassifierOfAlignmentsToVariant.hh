//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include <iostream>
#include <vector>

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

#include "common/CountTable.hh"

namespace ehunter {


class ClassifierOfAlignmentsToVariant
{
public:
    static const graphtools::NodeId kInvalidNodeId;

    ClassifierOfAlignmentsToVariant(std::vector<graphtools::NodeId> targetNodes);

    void classify(const graphtools::GraphAlignment& graphAlignment);

    const CountTable& countsOfReadsFlankingUpstream() const { return countsOfReadsFlankingUpstream_; }
    const CountTable& countsOfReadsFlankingDownstream() const { return countsOfReadsFlankingDownstream_; }
    const CountTable& countsOfSpanningReads() const { return countsOfSpanningReads_; }
    int numBypassingReads() const { return numBypassingReads_; }

private:
    std::vector<graphtools::NodeId> targetNodes_;
    graphtools::NodeId firstBundleNode_;
    graphtools::NodeId lastBundleNode_;

    CountTable countsOfReadsFlankingUpstream_;
    CountTable countsOfReadsFlankingDownstream_;
    CountTable countsOfSpanningReads_;
    int numBypassingReads_ = 0;
};

}
