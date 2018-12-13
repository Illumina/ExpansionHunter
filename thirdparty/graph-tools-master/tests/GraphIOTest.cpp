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

#include <iostream>
#include <mutex>

#include "nlohmann/json.hpp"
#include "gtest/gtest.h"
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "graphIO/BamWriter.hh"
#include "graphIO/GraphJson.hh"
#include "graphIO/ReferenceGenome.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"

#include <stdio.h>
#include <string.h>

using std::string;
namespace fs = boost::filesystem;
using namespace testing;
using namespace graphtools;
using namespace graphIO;

class ReferenceGenome : public testing::Test
{
protected:
    // Write Fasta file shared by all refGenome dependent tests.
    static void SetUpTestCase()
    {
        auto tmpFasta = fs::unique_path(fs::temp_directory_path() / "%%%%%%%%.fa");
        fs::ofstream fastaOut(tmpFasta);
        fastaOut << ">chr12\n"
                 << "AAAAAGGGGG" << std::endl;
        fastaOut.close();
        fastaPath = tmpFasta.string();
    }

    static void TearDownTestCase()
    {
        auto tmpFasta = fs::path(fastaPath);
        fs::remove(tmpFasta);
        fs::remove(tmpFasta.replace_extension(".fa.fai"));
    }

    static string fastaPath;
};
string ReferenceGenome::fastaPath = "";

TEST_F(ReferenceGenome, GetSequence_success)
{
    RefGenome ref(fastaPath);
    string const seq = ref.extractSeq(ReferenceInterval("chr12", 3, 6));
    ASSERT_EQ("AAG", seq);
}

TEST_F(ReferenceGenome, ParseInvalidRegion_throws)
{
    RefGenome ref(fastaPath);
    EXPECT_ANY_THROW(ReferenceInterval::parseRegion("chr12-4-6"));
}

TEST_F(ReferenceGenome, NonExistingSequence_throws)
{
    RefGenome ref(fastaPath);
    EXPECT_ANY_THROW(ref.extractSeq(ReferenceInterval("chr12", 4, 11)));
    EXPECT_ANY_THROW(ref.extractSeq(ReferenceInterval("chr13", 4, 6)));
}

TEST_F(ReferenceGenome, ParseRegion_success)
{
    RefGenome ref(fastaPath);
    ReferenceInterval reg = ReferenceInterval::parseRegion("chr12:4-6");

    ASSERT_EQ("chr12", reg.contig);
    ASSERT_EQ(4, reg.start);
    ASSERT_EQ(6, reg.end);
}

std::mutex HtsLibMutex;
class BamWriterTest : public Test
{
public:
    fs::path bamFile;
    void SetUp() override
    {
        HtsLibMutex.lock();
        bamFile = fs::unique_path(fs::temp_directory_path() / "%%%%%%%%.bam");
    }

    void TearDown() override
    {
        HtsLibMutex.unlock();
        fs::remove(bamFile);
    }
};

TEST_F(BamWriterTest, UnplacedAlignment_SingleRead)
{
    ReferenceContigs contigs = {};
    BamWriter bw(bamFile.string(), contigs);
    auto aln = bw.makeAlignment("Read2", "GATC", std::vector<int>(), BamWriter::PairingInfo::Unpaired, "1(1M2D)2(4M)");
    EXPECT_NO_THROW(bw.writeAlignment(aln));
    HtsLibMutex.unlock();
    ASSERT_EQ("", aln.chromName);
    ASSERT_EQ(-1, aln.pos);
    ASSERT_EQ(false, aln.isMate1);
    ASSERT_EQ(false, aln.isPaired);
    // TODO GT-538: Add test to validate BAM
}

TEST_F(BamWriterTest, UnplacedAlignment_PairedReads)
{
    ReferenceContigs contigs = {};
    BamWriter bw(bamFile.string(), contigs);
    auto aln1 = bw.makeAlignment("Read1", "ATTAC", std::vector<int>(), BamWriter::PairingInfo::FirstMate, "1(3M)");
    ASSERT_NO_THROW(bw.writeAlignment(aln1));
    ASSERT_TRUE(aln1.isMate1);
    ASSERT_TRUE(aln1.isPaired);
    auto aln2
        = bw.makeAlignment("Read1", "GATC", std::vector<int>(), BamWriter::PairingInfo::SecondMate, "1(1M2D)2(4M)");
    ASSERT_NO_THROW(bw.writeAlignment(aln2));
    ASSERT_FALSE(aln2.isMate1);
    ASSERT_TRUE(aln2.isPaired);
    // TODO GT-538: Add test to validate BAM
}

TEST_F(BamWriterTest, PlacedAlignment_SingleRead)
{
    ReferenceContigs contigs = { std::make_pair("chr1", 10), std::make_pair("chr2", 20) };
    BamWriter bw(bamFile.string(), contigs);
    Graph graph = makeSwapGraph("AAAA", "C", "T", "GGGG");
    GraphReferenceMapping mapping(&graph);
    mapping.addMapping(0, ReferenceInterval("chr2", 10, 14));
    Path path(&graph, 2, { 0, 1, 3 }, 3);
    std::vector<Alignment> alignments = { Alignment(2, "2M"), Alignment(0, "1M"), Alignment(0, "3M") };
    GraphAlignment gAlign(path, alignments);

    auto aln = bw.makeAlignment(mapping, "read1", "AACGGG", {}, BamWriter::PairingInfo::Unpaired, gAlign);
    ASSERT_NO_THROW(bw.writeAlignment(aln));
    ASSERT_EQ("chr2", aln.chromName);
    ASSERT_EQ(12, aln.pos);
    ASSERT_EQ("0[Ref start: 2, 2M]1[Ref start: 0, 1M]3[Ref start: 0, 3M]", aln.graphCigar);
}

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

TEST_F(ReferenceGenome, LoadGraphSequence_Success)
{
    Json jGraph;
    jGraph["reference_genome"] = fastaPath;
    jGraph["nodes"] = { { { "name", "n1" }, { "reference", "chr12:3-7" } } };
    jGraph["edges"] = Json::array();

    Graph const graph = parseGraph(jGraph);

    ASSERT_EQ("AAGG", graph.nodeSeq(0));
}

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