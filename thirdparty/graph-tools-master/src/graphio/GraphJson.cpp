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

#include "graphio/GraphJson.hh"

#include <boost/optional.hpp>
#include <fstream>

using std::string;

namespace graphtools
{

Graph loadGraph(string const& jsonPath)
{
    std::ifstream jsonFile(jsonPath);
    assert(jsonFile.good());
    Json json;
    jsonFile >> json;
    return parseGraph((json.count("graph") == 1) ? json["graph"] : json);
}

Graph parseGraph(Json const& jsonGraph)
{
    // boost::optional<RefGenome> genome;
    // const auto refFasta = jsonGraph.find("reference_genome");
    // if (refFasta != jsonGraph.end())
    //{
    //    genome.emplace(refFasta->get<string>());
    //}

    Json::array_t nodes = jsonGraph.value("nodes", Json::array());
    auto const nNodes = nodes.size();
    string const& graphId = jsonGraph.value("graph_id", "");
    Graph graph(nNodes, graphId);

    std::unordered_map<string, NodeId> nodeIds; // NodeName -> NameID
    NodeId nodeIndex = 0;
    for (auto const& jsonNode : nodes)
    {
        string const name = jsonNode.at("name");
        graph.setNodeName(nodeIndex, name);
        auto const seq = jsonNode.find("sequence");
        if (seq != jsonNode.end())
        {
            graph.setNodeSeq(nodeIndex, *seq);
        }
        else
        {
            throw std::runtime_error("Node has an invalid sequence: " + graph.nodeName(nodeIndex));
            // auto const refRegion = jsonNode.find("reference");
            // if (refRegion == jsonNode.end())
            //{
            //    throw std::runtime_error("Node has no sequence: " + graph.nodeName(nodeIndex));
            //}
            // if (!genome)
            //{
            //    throw std::runtime_error("Need 'referenceGenome' FASTA file to use reference nodes");
            //}
            // auto const interval = ReferenceInterval::parseRegion(*refRegion);
            // graph.setNodeSeq(nodeIndex, genome->extractSeq(interval));
        }
        assert(nodeIds.count(name) == 0);
        nodeIds[name] = nodeIndex;
        nodeIndex++;
    }
    for (auto const& jsonEdge : jsonGraph.value("edges", Json::array()))
    {
        NodeId const sourceNode = nodeIds.at(jsonEdge.at("from"));
        NodeId const sinkNode = nodeIds.at(jsonEdge.at("to"));
        graph.addEdge(sourceNode, sinkNode);
        const auto labels = jsonEdge.find("labels");
        if (labels != jsonEdge.end())
        {
            for (string const& label : *labels)
            {
                graph.addLabelToEdge(sourceNode, sinkNode, label);
            }
        }
    }
    return graph;
}

Json graphToJson(Graph const& graph)
{
    Json json;
    if (!graph.graphId.empty())
    {
        json["graph_id"] = graph.graphId;
    }
    json["nodes"] = Json::array();
    for (size_t i = 0; i != graph.numNodes(); ++i)
    {
        json["nodes"].push_back({ { "name", graph.nodeName(i) }, { "sequence", graph.nodeSeq(i) } });
    }
    json["edges"] = Json::array();
    for (size_t n1 = 0; n1 < graph.numNodes(); ++n1)
    {
        for (NodeId n2 : graph.successors(n1))
        {
            Json edge = { { "from", graph.nodeName(n1) }, { "to", graph.nodeName(n2) } };
            auto const& labels = graph.edgeLabels(n1, n2);
            if (labels.size() > 0)
            {
                edge["labels"] = labels;
            }

            json["edges"].push_back(std::move(edge));
        }
    }
    return json;
}

GraphReferenceMapping parseReferenceMapping(Json const& jRefmap, Graph const& graph)
{
    Json::array_t nodes = jRefmap.value("nodes", Json::array());

    GraphReferenceMapping refmap(&graph);
    NodeId nodeIndex = 0;
    for (auto const& jNode : nodes)
    {
        const auto refInterval = jNode.find("reference");
        if (refInterval != jNode.end())
        {
            auto const region = ReferenceInterval::parseRegion(*refInterval);
            refmap.addMapping(nodeIndex, region);
        }
        ++nodeIndex;
    }
    return refmap;
}
}