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

#include "graphcore/GraphBuilders.hh"

#include "gtest/gtest.h"

using std::string;

using namespace graphtools;

TEST(CreatingGraphs, TypicalSequences_DeletionGraphCreated)
{
    const string left_flank = "AATT";
    const string deletion = "CCCC";
    const string right_flank = "GGGCC";
    Graph graph = makeDeletionGraph(left_flank, deletion, right_flank);

    EXPECT_EQ(3ul, graph.numNodes());
    EXPECT_EQ(left_flank, graph.nodeSeq(0));
    EXPECT_EQ(deletion, graph.nodeSeq(1));
    EXPECT_EQ(right_flank, graph.nodeSeq(2));
    EXPECT_TRUE(graph.hasEdge(0, 1));
    EXPECT_TRUE(graph.hasEdge(0, 2));
    EXPECT_TRUE(graph.hasEdge(1, 2));
}

TEST(CreatingGraphs, TypicalSequences_SwapGraphCreated)
{
    const string left_flank = "AATT";
    const string deletion = "CCCC";
    const string insertion = "TTTT";
    const string right_flank = "GGGCC";
    Graph graph = makeSwapGraph(left_flank, deletion, insertion, right_flank);

    EXPECT_EQ(4ul, graph.numNodes());
    EXPECT_EQ(left_flank, graph.nodeSeq(0));
    EXPECT_EQ(deletion, graph.nodeSeq(1));
    EXPECT_EQ(insertion, graph.nodeSeq(2));
    EXPECT_EQ(right_flank, graph.nodeSeq(3));
    EXPECT_TRUE(graph.hasEdge(0, 1));
    EXPECT_TRUE(graph.hasEdge(0, 2));
    EXPECT_TRUE(graph.hasEdge(1, 3));
    EXPECT_TRUE(graph.hasEdge(2, 3));
}

TEST(CreatingGraphs, TypicalSequences_DoubleSwapGraphCreated)
{
    const string left_flank = "AATT";
    const string deletion1 = "CCCC";
    const string insertion1 = "TTTT";
    const string middle = "CCCC";
    const string deletion2 = "AAAA";
    const string insertion2 = "GGGG";
    const string right_flank = "GGGCC";
    Graph graph = makeDoubleSwapGraph(left_flank, deletion1, insertion1, middle, deletion2, insertion2, right_flank);

    EXPECT_EQ(7ul, graph.numNodes());
    EXPECT_EQ(left_flank, graph.nodeSeq(0));
    EXPECT_EQ(deletion1, graph.nodeSeq(1));
    EXPECT_EQ(insertion1, graph.nodeSeq(2));
    EXPECT_EQ(middle, graph.nodeSeq(3));
    EXPECT_EQ(deletion2, graph.nodeSeq(4));
    EXPECT_EQ(insertion2, graph.nodeSeq(5));
    EXPECT_EQ(right_flank, graph.nodeSeq(6));
    EXPECT_TRUE(graph.hasEdge(0, 1));
    EXPECT_TRUE(graph.hasEdge(0, 2));
    EXPECT_TRUE(graph.hasEdge(1, 3));
    EXPECT_TRUE(graph.hasEdge(2, 3));
    EXPECT_TRUE(graph.hasEdge(3, 4));
    EXPECT_TRUE(graph.hasEdge(3, 5));
    EXPECT_TRUE(graph.hasEdge(4, 6));
    EXPECT_TRUE(graph.hasEdge(5, 6));
}

TEST(CreatingGraphs, TypicalSequences_LooplessStrGraphCreated)
{
    const string left_flank = "AATT";
    const string repeat_unit = "CGG";
    const string right_flank = "ATTT";
    const int32_t read_len = 10;
    Graph graph = makeLooplessStrGraph(read_len, left_flank, repeat_unit, right_flank);

    ASSERT_EQ(6ul, graph.numNodes());
    EXPECT_EQ(left_flank, graph.nodeSeq(0));
    EXPECT_EQ(repeat_unit, graph.nodeSeq(1));
    EXPECT_EQ(repeat_unit, graph.nodeSeq(2));
    EXPECT_EQ(repeat_unit, graph.nodeSeq(3));
    EXPECT_EQ(repeat_unit, graph.nodeSeq(4));
    EXPECT_EQ(right_flank, graph.nodeSeq(5));

    EXPECT_TRUE(graph.hasEdge(0, 5));

    EXPECT_TRUE(graph.hasEdge(0, 1));
    EXPECT_TRUE(graph.hasEdge(1, 5));

    EXPECT_TRUE(graph.hasEdge(1, 2));
    EXPECT_TRUE(graph.hasEdge(2, 5));

    EXPECT_TRUE(graph.hasEdge(2, 3));
    EXPECT_TRUE(graph.hasEdge(3, 5));

    EXPECT_TRUE(graph.hasEdge(3, 4));
    EXPECT_TRUE(graph.hasEdge(4, 5));
}

TEST(CreatingGraphs, TypicalSequences_StrGraphCreated)
{
    const string left_flank = "AATT";
    const string repeat_unit = "CGG";
    const string right_flank = "ATTT";
    Graph graph = makeStrGraph(left_flank, repeat_unit, right_flank);

    ASSERT_EQ(3ul, graph.numNodes());
    EXPECT_EQ(left_flank, graph.nodeSeq(0));
    EXPECT_EQ(repeat_unit, graph.nodeSeq(1));
    EXPECT_EQ(right_flank, graph.nodeSeq(2));

    EXPECT_TRUE(graph.hasEdge(0, 1));
    EXPECT_TRUE(graph.hasEdge(0, 2));
    EXPECT_TRUE(graph.hasEdge(1, 1));
    EXPECT_TRUE(graph.hasEdge(1, 2));
}
