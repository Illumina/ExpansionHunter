//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Roman Petrovski <RPetrovski@illumina.com>
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <stdio.h>
#include <string.h>

#include "graphalign/DagAlignerAffine.hh"
#include "graphalign/GappedAligner.hh"
#include "graphalign/PinnedDagAligner.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphcore/Path.hh"

#include "gtest/gtest.h"

using std::string;
using namespace testing;
using namespace graphtools;
using namespace graphalign;
using namespace graphalign::dagAligner;

template <typename T> std::string toString(const T& obj)
{
    std::stringstream ss;
    ss << obj;
    return ss.str();
}

TEST(SimpleGraph, AlignSuffix)
{
    Graph graph = makeSwapGraph("AAAA", "C", "T", "GGGG");
    graph.addEdge(1, 1);
    Path seed(
        &graph, 1,
        {
            0,
        },
        3);
    PinnedDagAligner dag_pinned_aligner(1, -1, 0, -2);
    int32_t top_dag_score = INT32_MIN;
    std::list<std::pair<Path, Alignment>> res = dag_pinned_aligner.prefixAlign(seed, "ACGG", 8, top_dag_score);
    EXPECT_EQ(4, top_dag_score);
    EXPECT_EQ("(0@1)-(1)-(3@2)", toString(res.front().first));
    EXPECT_EQ("4M", res.front().second.generateCigar());
}

TEST(SimpleGraph, AlignPrefix)
{
    Graph graph = makeSwapGraph("AAAA", "C", "T", "GGGG");
    graph.addEdge(1, 1);
    Path seed(&graph, 0, { 2, 3 }, 2);
    PinnedDagAligner dag_pinned_aligner(1, -1, 0, -2);
    int32_t top_dag_score = INT32_MIN;
    std::list<std::pair<Path, Alignment>> res = dag_pinned_aligner.suffixAlign(seed, "AA", 8, top_dag_score);
    EXPECT_EQ(2, top_dag_score);
    EXPECT_EQ("(0@2)-(2)-(3@2)", toString(res.front().first));
    EXPECT_EQ("2M", res.front().second.generateCigar());
}
