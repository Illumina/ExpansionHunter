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
