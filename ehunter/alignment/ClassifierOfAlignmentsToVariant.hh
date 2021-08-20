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
#include <vector>

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

#include "core/CountTable.hh"

namespace ehunter
{

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
