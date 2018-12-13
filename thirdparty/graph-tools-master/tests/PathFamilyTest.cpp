//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>
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

#include "gtest/gtest.h"

#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphcore/Path.hh"
#include "graphcore/PathFamily.hh"
#include <iostream>

using std::string;
using std::vector;

using namespace graphtools;

TEST(Creation, AddingEdges_ExpectedSize)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    PathFamily family(&graph);
    ASSERT_EQ(0u, family.edges().size());
    family.addEdge(0, 1);
    ASSERT_EQ(1u, family.edges().size());
}

TEST(Creation, FromLabel_Edgeset)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    graph.addLabelToEdge(0, 2, "foo");
    PathFamily family(&graph, "foo");
    ASSERT_EQ(1u, family.edges().size());
    family.addEdge(1, 2);
    family.setLabel("foo");
    PathFamily family2(&graph, "foo");
    ASSERT_EQ(family, family2);
}

TEST(Creation, CopyConstructor_Independent)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    PathFamily fam1(&graph);
    fam1.addEdge(0, 1);
    PathFamily fam2(fam1);
    ASSERT_EQ(fam1, fam2);
    fam1.addEdge(0, 2);
    ASSERT_NE(fam1, fam2);
}

TEST(Paths, ContainsPath_MatchFull)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    PathFamily family(&graph);
    family.addEdge(0, 1);
    family.addEdge(1, 2);
    Path path(&graph, 0, { 0, 1, 2 }, 0);
    ASSERT_TRUE(family.containsPath(path));
    Path path2(&graph, 0, { 0, 2 }, 0);
    ASSERT_FALSE(family.containsPath(path2));
}

TEST(Paths, ContainsPath_MatchPartialPath)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    PathFamily family(&graph);
    family.addEdge(0, 1);
    family.addEdge(1, 2);
    Path path(&graph, 0, { 0, 1 }, 0);
    ASSERT_TRUE(family.containsPath(path));
    Path path2(&graph, 0, { 1, 2 }, 0);
    ASSERT_TRUE(family.containsPath(path2));
}

TEST(Paths, ContainsPath_MatchPartialFamily)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    PathFamily family(&graph);
    family.addEdge(0, 1);
    Path path(&graph, 0, { 0, 1, 2 }, 0);
    ASSERT_TRUE(family.containsPath(path));
    Path path2(&graph, 0, { 1, 2 }, 0);
    ASSERT_FALSE(family.containsPath(path2));
}

TEST(Paths, ContainsPath_MatchAmbiguous)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    PathFamily family(&graph);
    family.addEdge(0, 1);
    family.addEdge(0, 2);
    Path path2(&graph, 0, { 0, 2 }, 0);
    ASSERT_TRUE(family.containsPath(path2));
    Path path(&graph, 0, { 0, 1, 2 }, 0);
    ASSERT_FALSE(family.containsPath(path));
}

TEST(Compare, AddingEdges_Equality)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    PathFamily family(&graph);
    family.addEdge(0, 1);
    PathFamily family2(&graph);
    family2.addEdge(0, 1);
    ASSERT_EQ(family, family2);
    family.addEdge(0, 2);
    ASSERT_NE(family, family2);
}
TEST(Compare, AddingEdges_Includes)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    graph.addLabelToEdge(0, 1, "foo");
    PathFamily fam1(&graph, "foo");
    PathFamily fam2(&graph, "foo");
    ASSERT_TRUE(fam1.includes(fam2));
    fam1.addEdge(0, 2);
    ASSERT_TRUE(fam1.includes(fam2));
    ASSERT_FALSE(fam2.includes(fam1));
}