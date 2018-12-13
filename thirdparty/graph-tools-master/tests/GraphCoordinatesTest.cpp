//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Peter Krusche <pkrusche@illumina.com>
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

#include "graphcore/GraphCoordinates.hh"

#include "gtest/gtest.h"

using namespace testing;
using namespace graphtools;

class GraphCoordinatesTest : public Test
{
public:
    Graph graph{ 4 };

    void SetUp() override
    {
        graph.setNodeName(0, "LF");
        graph.setNodeSeq(0, "AAAAAAAAAAA");

        graph.setNodeName(1, "P1");
        graph.setNodeSeq(1, "TTTTTT");

        graph.setNodeName(2, "Q1");
        graph.setNodeSeq(2, "GGGGGGGG");

        graph.setNodeName(3, "RF");
        graph.setNodeSeq(3, "AAAAAAAAAAA");

        /**
         *
         * LF           RF
         *  |           ^
         *  |           |
         *  *-> P1 -----*
         *  |           |
         *  *-> Q1 -----*
         */

        graph.addEdge(0, 1);
        graph.addEdge(0, 2);
        graph.addEdge(1, 3);
        graph.addEdge(2, 3);
    }
};

TEST_F(GraphCoordinatesTest, CanonicalPositionLookup)
{
    GraphCoordinates coordinates(&graph);

    // LF has offset 0
    ASSERT_EQ(static_cast<uint64_t>(6), coordinates.canonicalPos("LF", 6));
    ASSERT_EQ(static_cast<uint64_t>(11 + 4), coordinates.canonicalPos("P1", 4));
    ASSERT_EQ(static_cast<uint64_t>(11 + 6 + 3), coordinates.canonicalPos("Q1", 3));
    ASSERT_EQ(static_cast<uint64_t>(11 + 6 + 8 + 2), coordinates.canonicalPos("RF", 2));
}

TEST_F(GraphCoordinatesTest, ReverseLookup)
{
    GraphCoordinates coordinates(&graph);
    for (size_t j = 0; j < graph.nodeSeq(0).size(); ++j) // node 0 == LF
    {
        std::string n;
        uint64_t offset = static_cast<uint64_t>(-1);
        coordinates.nodeAndOffset(j, n, offset);
        ASSERT_EQ("LF", n);
        ASSERT_EQ(static_cast<uint64_t>(j), offset);
    }

    for (size_t j = 0; j < graph.nodeSeq(1).size(); ++j) // node 1 == P1
    {
        std::string n;
        uint64_t offset = static_cast<uint64_t>(-1);
        coordinates.nodeAndOffset(11 + j, n, offset);
        ASSERT_EQ("P1", n);
        ASSERT_EQ(static_cast<uint64_t>(j), offset);
    }
    for (size_t j = 0; j < graph.nodeSeq(2).size(); ++j) // node 2 == Q1
    {
        std::string n;
        uint64_t offset = static_cast<uint64_t>(-1);
        coordinates.nodeAndOffset(11 + 6 + j, n, offset);
        ASSERT_EQ("Q1", n);
        ASSERT_EQ(static_cast<uint64_t>(j), offset);
    }
    for (size_t j = 0; j < graph.nodeSeq(3).size(); ++j) // node 3 == RF
    {
        std::string n;
        uint64_t offset = static_cast<uint64_t>(-1);
        coordinates.nodeAndOffset(11 + 6 + 8 + j, n, offset);
        ASSERT_EQ("RF", n);
        ASSERT_EQ(static_cast<uint64_t>(j), offset);
    }
}

TEST_F(GraphCoordinatesTest, DistanceComputation)
{
    GraphCoordinates coordinates(&graph);

    // both on LF
    ASSERT_EQ(static_cast<uint64_t>(5), coordinates.distance(10, 5));
    ASSERT_EQ(static_cast<uint64_t>(5), coordinates.distance(5, 10));

    // one on LF, one on neighbour (P1 or Q1)
    ASSERT_EQ(static_cast<uint64_t>(8), coordinates.distance(14, 6));
    ASSERT_EQ(static_cast<uint64_t>(8), coordinates.distance(20, 6));

    // LF -> RF should go via P1 because this is shorter
    ASSERT_EQ(static_cast<uint64_t>(9 + 6 + 4), coordinates.distance(2, 11 + 6 + 8 + 4));
}
