//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>
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
