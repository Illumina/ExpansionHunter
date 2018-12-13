//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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
