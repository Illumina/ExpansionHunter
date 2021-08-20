//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>
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

#include "graphcore/GraphReferenceMapping.hh"
#include "graphcore/GraphBuilders.hh"

#include <memory>

#include "gtest/gtest.h"

#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

using namespace graphtools;

class GraphMapping : public testing::Test
{
public:
    GraphMapping()
        : graph(makeSwapGraph("AAAA", "C", "T", "GGGG"))
        , mapping(&graph)
    {
        mapping.addMapping(0, ReferenceInterval("chr1", 10, 14));
        mapping.addMapping(2, ReferenceInterval("chr1", 15, 16));
        mapping.addMapping(3, ReferenceInterval("chr1", 16, 20));
    }

    Graph graph;
    GraphReferenceMapping mapping;
};
TEST_F(GraphMapping, MapNodePosition_Success)
{
    ASSERT_EQ(ReferenceInterval::makePosition("chr1", 10), mapping.map(0, 0));
    ASSERT_EQ(ReferenceInterval::makePosition("chr1", 13), mapping.map(0, 3));
    ASSERT_EQ(ReferenceInterval::makePosition("chr1", 15), mapping.map(2, 0));
}
TEST_F(GraphMapping, UnmappedNode_ReturnEmpty) { ASSERT_FALSE(mapping.map(1, 0)); }
TEST_F(GraphMapping, MapInvalidPos_Throws)
{
    ASSERT_ANY_THROW(mapping.map(2, 1)); // Outside Node
    ASSERT_ANY_THROW(mapping.map(5, 0)); // Invalid node
}
TEST_F(GraphMapping, MapPath_StartingNode)
{
    auto const map = mapping.map(Path(&graph, 1, { 0, 1, 3 }, 4));
    ASSERT_EQ(ReferenceInterval::makePosition("chr1", 11), map);
}
TEST_F(GraphMapping, MapPath_ExtendingNode)
{
    auto const map = mapping.map(Path(&graph, 0, { 1, 3 }, 4));
    ASSERT_EQ(ReferenceInterval::makePosition("chr1", 16), map);
}
TEST_F(GraphMapping, UnmappedPath_ReturnEmpty)
{
    auto const map = mapping.map(Path(&graph, 0, { 1 }, 1));
    ASSERT_FALSE(map);
}
