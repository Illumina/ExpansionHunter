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

#include <iostream>

#include "nlohmann/json.hpp"
#include "gtest/gtest.h"

#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphio/GraphJson.hh"

using std::string;

TEST(GraphLoading, ValidGraph_Loaded)
{
    Json jGraph;
    jGraph["nodes"] = { { { "name", "n1" }, { "sequence", "AATG" } },
                        { { "name", "n2" }, { "sequence", "AA" } },
                        { { "name", "n3" }, { "sequence", "TG" } } };
    jGraph["edges"] = { { { "from", "n1" }, { "to", "n2" } },
                        { { "from", "n2" }, { "to", "n3" } },
                        { { "from", "n2" }, { "to", "n2" } } };
    jGraph["graph_id"] = "TestGraph";

    Graph const graph = parseGraph(jGraph);

    ASSERT_EQ("TestGraph", graph.graphId);
    ASSERT_EQ(jGraph["nodes"].size(), graph.numNodes());
    ASSERT_EQ(jGraph["edges"].size(), graph.numEdges());
    for (size_t i = 0; i != graph.numNodes(); ++i)
    {
        ASSERT_EQ(jGraph["nodes"][i]["name"], graph.nodeName(i));
        ASSERT_EQ(jGraph["nodes"][i]["sequence"], graph.nodeSeq(i));
    }
    ASSERT_TRUE(graph.hasEdge(0, 1));
    ASSERT_TRUE(graph.hasEdge(1, 2));
    ASSERT_TRUE(graph.hasEdge(1, 1));
    ASSERT_FALSE(graph.hasEdge(0, 0));
    ASSERT_FALSE(graph.hasEdge(0, 2));
}

TEST(GraphLoading, MissingSequence_Throws)
{
    Json jGraph;
    jGraph["nodes"] = {
        { { "name", "n1" } },
    };
    jGraph["edges"] = Json::array();

    ASSERT_ANY_THROW(parseGraph(jGraph));
}

TEST(GraphLoading, EmptySequence_Throws)
{
    Json jGraph;
    jGraph["nodes"] = {
        { { "name", "n1" }, { "sequence", "" } },
    };
    jGraph["edges"] = Json::array();

    ASSERT_ANY_THROW(parseGraph(jGraph));
}

TEST(GraphLoading, InvalidEdgeNode_Throws)
{
    Json jGraph;
    jGraph["nodes"] = {
        { { "name", "n1" }, { "sequence", "AATG" } },
    };
    jGraph["edges"] = {
        { { "from", "n1" }, { "to", "n2" } },
    };

    ASSERT_ANY_THROW(parseGraph(jGraph));
}

TEST(GraphLoading, BackwardsEdge_Throws)
{
    Json jGraph;
    jGraph["nodes"] = {
        { { "name", "n1" }, { "sequence", "AATG" } },
        { { "name", "n2" }, { "sequence", "AATG" } },
    };
    jGraph["edges"] = {
        { { "from", "n2" }, { "to", "n1" } },
    };

    ASSERT_ANY_THROW(parseGraph(jGraph));
}

/*

TEST(ReferenceGenome, LoadGraphSequence_Success)
{
    Json jGraph;
    jGraph["reference_genome"] = fastaPath;
    jGraph["nodes"] = { { { "name", "n1" }, { "reference", "chr12:3-7" } } };
    jGraph["edges"] = Json::array();

    Graph const graph = parseGraph(jGraph);

    ASSERT_EQ("AAGG", graph.nodeSeq(0));
}

*/

TEST(GraphLoading, MissingReference_Throws)
{
    Json jGraph;
    jGraph["nodes"] = { { { "name", "n1" }, { "reference", "chr12:4-7" } } };
    jGraph["edges"] = Json::array();

    ASSERT_ANY_THROW(parseGraph(jGraph));
}

TEST(GraphWriting, EmptyGraph_RoundTrip)
{
    Graph graph(0);
    Json jGraph = graphToJson(graph);
    Graph newGraph = parseGraph(jGraph);

    ASSERT_EQ((size_t)0, newGraph.numNodes());
}

TEST(GraphWriting, Graph_RoundTrip)
{
    Graph graph(2, "Small Graph");
    graph.setNodeName(0, "n0");
    graph.setNodeSeq(0, "AA");
    graph.setNodeName(1, "n1");
    graph.setNodeSeq(1, "TT");
    graph.addEdge(0, 1);
    graph.addEdge(1, 1);
    graph.addLabelToEdge(1, 1, "foo");

    Json jGraph = graphToJson(graph);
    Graph newGraph = parseGraph(jGraph);

    ASSERT_EQ("Small Graph", graph.graphId);
    ASSERT_EQ(graph.numNodes(), newGraph.numNodes());
    ASSERT_EQ(graph.numEdges(), newGraph.numEdges());
    for (size_t i = 0; i != graph.numNodes(); ++i)
    {
        ASSERT_EQ(graph.nodeName(i), newGraph.nodeName(i));
        ASSERT_EQ(graph.nodeSeq(i), newGraph.nodeSeq(i));
    }
    ASSERT_TRUE(newGraph.hasEdge(0, 1));
    ASSERT_TRUE(newGraph.hasEdge(1, 1));
    ASSERT_FALSE(newGraph.hasEdge(0, 0));
    ASSERT_EQ(graph.edgeLabels(1, 1), newGraph.edgeLabels(1, 1));
}

/*

TEST_F(ReferenceGenome, LoadGraphMapping_Success)
{
    Json jGraph;
    jGraph["reference_genome"] = fastaPath;
    jGraph["nodes"]
        = { { { "name", "n1" }, { "reference", "chr12:4-7" } }, { { "name", "n2" }, { "sequence", "TCGA" } } };
    jGraph["edges"] = Json::array();

    Graph const graph = parseGraph(jGraph);
    GraphReferenceMapping const refmap = parseReferenceMapping(jGraph, graph);

    auto const pos = refmap.map(0, 2);
    ASSERT_TRUE(pos);
    ASSERT_EQ("chr12", pos->contig);
    ASSERT_EQ(6, pos->start);

    ASSERT_FALSE(refmap.map(1, 2)); // Node without mapping
    ASSERT_ANY_THROW(refmap.map(0, 3)); // Position outside node
    ASSERT_ANY_THROW(refmap.map(2, 0)); // Nonexistent node
}

*/
