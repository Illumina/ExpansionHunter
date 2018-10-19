// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.

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
