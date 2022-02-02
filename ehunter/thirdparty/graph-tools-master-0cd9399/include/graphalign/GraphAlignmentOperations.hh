//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#pragma once

#include <cstdint>
#include <string>

#include "graphalign/GraphAlignment.hh"
#include "graphalign/LinearAlignment.hh"
#include "graphcore/Graph.hh"

namespace graphtools
{

// Checks if graph alignment is consistent with the given query sequence
bool checkConsistency(const GraphAlignment& graph_alignment, const std::string& query);

// Returns true if (disregarding soft clipping) the alignment starts and ends with a match
bool isLocalAlignment(const GraphAlignment& graph_alignment);

// A node CIGAR is a string of the form "<node_id>[linear alignment CIGAR]". This function extracts node_id and
// linear alignment CIGAR from the node CIGAR.
void splitNodeCigar(const std::string& node_cigar, std::string& cigar, NodeId& node_id);

/**
 * Converts graph CIGAR string to a graph alignment
 *
 * @param first_node_start Start position of the alignment on the first node
 * @param graph_cigar Graph CIGAR string
 * @param query Query sequence
 * @param graph_ptr Pointer to the graph
 * @return GraphAlignment
 */
GraphAlignment decodeGraphAlignment(int32_t first_node_start, const std::string& graph_cigar, const Graph* graph_ptr);

/**
 * Convert linear alignment to graph alignment by projecting it onto a compatible path
 *
 * @param linear_alignment: Linear alignment
 * @param path: Path to project the alignment onto; sequence of the path must be equal to alignment's reference sequence
 * @return Graph alignment composed of the input linear alignment and the (appropriately shrunk) path
 */
GraphAlignment projectAlignmentOntoGraph(Alignment linear_alignment, Path path);

// Splits query into pieces corresponding to each node that the alignment spans
std::list<std::string> getQuerySequencesForEachNode(const GraphAlignment& graph_alignment, const std::string& query);

// Encodes alignment as a three-row strings where the top corresponds to the query sequence, the bottom to the sequence
// of alignment's path, and the middle contains a "|" for each pair of matching bases; gaps are indicated by "-"; ends
// of nodes are denoted by ":"
std::string prettyPrint(const GraphAlignment& graph_alignment, const std::string& query);
}
