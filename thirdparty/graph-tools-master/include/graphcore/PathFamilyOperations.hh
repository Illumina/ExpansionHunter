//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "Path.hh"
#include "PathFamily.hh"

#include <limits>
#include <map>
#include <string>

namespace graphtools
{

/**
 * Generate path segments in a family which can be combined into longer paths.
 *
 * These segments are built by concatenating family edges into paths whenever
 * this is possible uniquely, ignoring repeat / self-loop edges.
 *
 * @param family the path family
 * @return a list of path segments
 */
std::list<Path> getPathSegmentsForFamily(graphtools::PathFamily const& family);

/**
 * Enumerate path segment combinations in family.
 *
 * Two path segments can be combined if they overlap or are adjacent and
 * (if adjacent on different nodes) their linking edge is in the family.
 *
 * @param family a path family
 * @param segments a set of path segments
 * @param[out] paths output list of paths
 * @param maxPaths maximum number of paths to generate
 * @return true if all paths were generated, false if maxPaths was used to limit the number of paths
 */
bool enumeratePathCombinationsInFamily(
    PathFamily const& family, std::list<Path> const& segments, std::list<Path>* paths, size_t maxPaths);

/**
 * Get all maximal paths in a path family, exhaustively
 *
 * Note that this function can generate an number of paths that is
 * exponential in the number of nodes.
 *
 * @param family the path family
 * @param[out] paths output list of paths
 * @param maxPaths limit the number of paths
 * @return true if all paths were enumerated, false if additional paths can be generated
 *
 * Note this will ignore self-edges / loops.
 *
 */
bool getMaximalPathsForFamily(
    graphtools::PathFamily const& family, std::list<Path>* paths, size_t maxPaths = std::numeric_limits<size_t>::max());

/**
 * Convert path to path family
 * @param graph must match the graph used by path
 * @param path a path on graph
 * @return a path family containing all edges in path
 */
graphtools::PathFamily pathToPathFamily(graphtools::Graph& graph, graphtools::Path const& path);

/**
 * Get all path families from edge labels on a graph
 * @param graph a graph (not const because path families are constructed from a non-const graph)
 * @return mapping of label -> path family
 */
std::map<std::string, graphtools::PathFamily> getPathFamiliesFromGraph(graphtools::Graph& graph);
}
