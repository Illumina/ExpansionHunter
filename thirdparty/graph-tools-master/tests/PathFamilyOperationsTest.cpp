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
#include "graphcore/PathFamilyOperations.hh"
#include <iostream>

using std::list;
using std::string;

using namespace graphtools;

TEST(PathsForPathFamily, GeneratePathsForPathFamily_DisjointPaths)
{
    Graph graph = makeDoubleSwapGraph("AAA", "CCCC", "GGG", "AAAA", "TTTT", "GG", "AA");
    PathFamily family(&graph);
    family.addEdge(0, 1);
    family.addEdge(1, 3);
    family.addEdge(5, 6);

    const list<Path> expected_paths{
        Path{ &graph, 0, { 0, 1, 3 }, 4 },
        Path{ &graph, 0, { 5, 6 }, 2 },
    };
    list<Path> observed_paths;
    ASSERT_TRUE(getMaximalPathsForFamily(family, &observed_paths));
    for (const auto& path : observed_paths)
    {
        ASSERT_TRUE(family.containsPath(path));
    }
    ASSERT_EQ(expected_paths, observed_paths);
}

TEST(PathsForPathFamily, GeneratePathsForPathFamily_LongPath)
{
    Graph graph = makeDoubleSwapGraph("AAA", "CCCC", "GGG", "AAAA", "TTTT", "GG", "AA");
    PathFamily family(&graph);
    family.addEdge(0, 1);
    family.addEdge(1, 3);
    family.addEdge(3, 4);
    family.addEdge(4, 6);

    const list<Path> expected_paths{
        Path{ &graph, 0, { 0, 1, 3, 4, 6 }, 2 },
    };
    list<Path> observed_paths;
    ASSERT_TRUE(getMaximalPathsForFamily(family, &observed_paths));
    for (const auto& path : observed_paths)
    {
        ASSERT_TRUE(family.containsPath(path));
    }
    ASSERT_EQ(expected_paths, observed_paths);
}

TEST(PathsForPathFamily, GeneratePathsForPathFamily_MultipleExtensions)
{
    Graph graph = makeDoubleSwapGraph("AAA", "CCCC", "GGG", "AAAA", "TTTT", "GG", "AA");
    PathFamily family(&graph);
    family.addEdge(1, 3);
    family.addEdge(2, 3);
    family.addEdge(3, 4);
    family.addEdge(3, 5);
    family.addEdge(4, 6);
    family.addEdge(5, 6);

    const list<Path> expected_paths{
        Path{ &graph, 0, { 1, 3, 4, 6 }, 2 },
        Path{ &graph, 0, { 1, 3, 5, 6 }, 2 },
        Path{ &graph, 0, { 2, 3, 4, 6 }, 2 },
        Path{ &graph, 0, { 2, 3, 5, 6 }, 2 },
    };
    list<Path> observed_paths;
    ASSERT_TRUE(getMaximalPathsForFamily(family, &observed_paths));

    for (const auto& path : observed_paths)
    {
        ASSERT_TRUE(family.containsPath(path));
    }

    ASSERT_EQ(expected_paths, observed_paths);
}

TEST(PathsForPathFamily, GeneratePathsForPathFamily_MultipleExtensionsSingleEdge)
{
    Graph graph(8);

    /*
       A      E
        \   /
         C=D
       /    \
      B      F
    */
    graph.setNodeName(0, "source");
    graph.setNodeName(1, "A");
    graph.setNodeName(2, "B");
    graph.setNodeName(3, "C");
    graph.setNodeName(4, "D");
    graph.setNodeName(5, "E");
    graph.setNodeName(6, "F");
    graph.setNodeName(7, "sink");

    graph.setNodeSeq(0, "N");
    graph.setNodeSeq(1, "A");
    graph.setNodeSeq(2, "A");
    graph.setNodeSeq(3, "A");
    graph.setNodeSeq(4, "A");
    graph.setNodeSeq(5, "A");
    graph.setNodeSeq(6, "A");
    graph.setNodeSeq(7, "N");

    graph.addEdge(0, 1);
    graph.addEdge(0, 2);
    graph.addEdge(1, 3);
    graph.addEdge(2, 3);
    graph.addEdge(3, 4);
    graph.addEdge(4, 5);
    graph.addEdge(4, 6);
    graph.addEdge(5, 7);
    graph.addEdge(6, 7);

    PathFamily family(&graph);
    family.addEdge(3, 4);

    const list<Path> expected_paths{
        Path{ &graph, 0, { 3, 4 }, 1 },
    };
    list<Path> observed_paths;
    ASSERT_TRUE(getMaximalPathsForFamily(family, &observed_paths));

    for (const auto& path : observed_paths)
    {
        ASSERT_TRUE(family.containsPath(path));
    }

    ASSERT_EQ(expected_paths, observed_paths);
}

TEST(PathsForPathFamily, GeneratePathsForPathFamily_LoopGraph)
{
    Graph graph = makeStrGraph("AAA", "TG", "CCC");
    PathFamily family(&graph);
    family.addEdge(0, 1);
    family.addEdge(1, 1);
    family.addEdge(1, 2);

    const list<Path> expected_paths{
        Path{ &graph, 0, { 0, 1, 2 }, 3 },
    };
    list<Path> observed_paths;
    ASSERT_TRUE(getMaximalPathsForFamily(family, &observed_paths, 5));
    for (const auto& path : observed_paths)
    {
        ASSERT_TRUE(family.containsPath(path));
    }
    ASSERT_EQ(expected_paths, observed_paths);
}

TEST(PathFamilyFromPath, GeneratePathFamilyFromPath_SimplePath)
{
    Graph graph = makeDoubleSwapGraph("A", "C", "T", "A", "G", "C", "T");
    const Path path{ &graph, 0, { 1, 3, 4 }, 0 };
    const PathFamily family = pathToPathFamily(graph, path);

    ASSERT_EQ(2ull, family.edges().size());
    ASSERT_NE(0u, family.edges().count({ 1, 3 }));
    ASSERT_NE(0u, family.edges().count({ 3, 4 }));
    ASSERT_TRUE(family.containsPath(path));
}

TEST(PathFamilyFromGraph, GeneratePathFamilyFromGraph_SimpleGraph)
{
    Graph graph = makeDoubleSwapGraph("A", "C", "T", "A", "G", "T", "C");
    graph.addLabelToEdge(0, 1, "A");
    graph.addLabelToEdge(1, 3, "A");
    graph.addLabelToEdge(3, 5, "A");
    graph.addLabelToEdge(5, 6, "A");
    graph.addLabelToEdge(0, 2, "B");
    graph.addLabelToEdge(2, 3, "B");
    graph.addLabelToEdge(3, 4, "B");
    graph.addLabelToEdge(4, 6, "B");

    const std::map<std::string, PathFamily> families = getPathFamiliesFromGraph(graph);

    ASSERT_EQ(2ull, families.size());
    ASSERT_NE(0u, families.count("A"));
    ASSERT_NE(0u, families.count("B"));

    const auto a_it = families.find("A");
    const auto b_it = families.find("B");
    ASSERT_EQ(4ull, a_it->second.edges().size());
    ASSERT_EQ(4ull, b_it->second.edges().size());

    ASSERT_NE(0u, a_it->second.edges().count({ 0, 1 }));
    ASSERT_NE(0u, a_it->second.edges().count({ 1, 3 }));
    ASSERT_NE(0u, a_it->second.edges().count({ 3, 5 }));
    ASSERT_NE(0u, a_it->second.edges().count({ 5, 6 }));

    ASSERT_NE(0u, b_it->second.edges().count({ 0, 2 }));
    ASSERT_NE(0u, b_it->second.edges().count({ 2, 3 }));
    ASSERT_NE(0u, b_it->second.edges().count({ 3, 4 }));
    ASSERT_NE(0u, b_it->second.edges().count({ 4, 6 }));
}

TEST(PathFamilyFromGraph, GeneratePathFamilyFromGraph_LoopGraph)
{
    Graph graph = makeStrGraph("A", "CT", "G");
    graph.addLabelToEdge(0, 2, "A");
    graph.addLabelToEdge(0, 1, "B");
    graph.addLabelToEdge(1, 1, "B");
    graph.addLabelToEdge(1, 2, "B");

    const std::map<std::string, PathFamily> families = getPathFamiliesFromGraph(graph);

    ASSERT_EQ(2ull, families.size());
    ASSERT_NE(0u, families.count("A"));
    ASSERT_NE(0u, families.count("B"));

    const auto a_it = families.find("A");
    const auto b_it = families.find("B");
    ASSERT_EQ(1ull, a_it->second.edges().size());
    ASSERT_EQ(3ull, b_it->second.edges().size());

    ASSERT_NE(0u, a_it->second.edges().count({ 0, 2 }));

    ASSERT_NE(0u, b_it->second.edges().count({ 0, 1 }));
    ASSERT_NE(0u, b_it->second.edges().count({ 1, 1 }));
    ASSERT_NE(0u, b_it->second.edges().count({ 1, 2 }));
}
