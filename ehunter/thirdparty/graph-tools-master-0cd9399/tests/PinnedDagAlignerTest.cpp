//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Roman Petrovski <RPetrovski@illumina.com>
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
