//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>
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

#include <list>

#include "graphcore/Path.hh"

namespace graphtools
{

/**
 * Computes all possible extensions of the path by the specified length in both directions
 *
 * @param start_extension_len Extension length of the start posision
 * @param end_extension_len Extension length of the end position
 * @return List of extended paths
 */
std::list<Path> extendPath(const Path& path, int32_t start_extension_len, int32_t end_extension_len);

// Computes all possible extensions of the start of the path by the specified length
std::list<Path> extendPathStart(const Path& path, int32_t extension_len);

// Computes all possible extensions of the end of the path by the specified length
std::list<Path> extendPathEnd(const Path& path, int32_t extension_len);

/**
 * Extend a path matching a query sequence to produce maximum exact + unique matches
 *
 * The path is extended at the end while the query string matches uniquely.
 *
 * @param path the path. We must have path.length() <= query.size() - qpos
 * @param query query string.
 * @param qpos position where path starts in query
 * @return extended path
 */
Path extendPathEndMatching(Path path, const std::string& query, size_t qpos);

/**
 * Extend a path matching a query sequence to produce maximum exact + unique matches
 *
 * The path is extended at the start while the query string matches uniquely.
 *
 * @param path the path. We must have path.length() <= query.size() - qpos
 * @param query query string.
 * @param[out] qpos position where path starts in query -- updated as path is extended at start
 * @return extended path
 */
Path extendPathStartMatching(Path path, const std::string& query, size_t& pos_in_query);

/**
 * Extend a path matching a query sequence to produce maximum exact + unique matches
 *
 * The path is extended at the start and end while the query string matches uniquely.
 *
 * @param path the path. We must have path.length() <= query.size() - qpos
 * @param query query string.
 * @param[out] qpos position where path starts in query -- updated as path is extended at start
 * @return extended path
 */
Path extendPathMatching(Path path, const std::string& query, size_t& pos_in_query);

/**
 * Splits sequence into segments corresponding to the path
 *
 * @param path Any path
 * @param sequence A string having the same length as the path
 * @return Segments of the sequence corresponding to nodes spanned by the path
 */
std::vector<std::string> splitSequenceByPath(const Path& path, const std::string& sequence);

/**
 * Return true if two paths are exactly adjacent
 * (i.e. p1 starts just before p2, or the other way around)
 *
 * @param p1 first path
 * @param p2 second path
 * @return true if the paths are adjacent
 */
bool checkIfPathsAdjacent(Path const& p1, Path const& p2);

/**
 * Return true if two paths overlap either prefix - suffix or suffix-prefix
 * @param p1 first path
 * @param p2 second path
 * @return true if the paths overlap
 */
bool checkPathPrefixSuffixOverlap(Path const& p1, Path const& p2);

/**
 * Paths can be merged if they overlap prefix-suffix / suffix-prefix.
 *
 * @param p1 first path
 * @param p2 second path
 * @return merged path
 */
Path mergePaths(Path const& p1, Path const& p2);

/**
 * Merge a set of paths
 *
 * This will merge paths iteratively until none of the resulting paths overlap
 *
 * @param paths a list of paths
 */
void greedyMerge(std::list<Path>& paths);

/**
 * Merge a set of paths
 *
 * This will merge paths exhaustively, each path is merged with all
 * paths it overlaps until we cannot merge anymore
 *
 * @param paths a list of paths
 */
void exhaustiveMerge(std::list<Path>& paths);

/**
 * Return the intersection(s) between two paths
 *
 * Multiple paths may be returned when paths diverge and re-join.
 *
 * @param p1 first path
 * @param p2 second path
 * @return merged path
 */
std::list<Path> intersectPaths(Path const& p1, Path const& p2);

/**
 * Returns subpaths corresponding to overlap of the input path with each node it passes through
 *
 * @param path: any path
 * @return single-node subpaths
 */
std::list<Path> generateSubpathForEachNode(const Path& path);

/**
 * Checks if two paths are bookended
 *
 * Two paths are considered bookended if positions of first path end and second path start are adjacent on the graph
 *
 * @param first_alignment: Any path
 * @param second_alignment: Any path on same graph
 * @return true if paths are bookended
 */
bool checkIfBookended(const Path& first_path, const Path& second_path);

/**
 * Concatenates two bookended paths into a longer path
 *
 * @param first_alignment: Any path
 * @param second_alignment: Any path that is bookended with the first
 * @return Concatenated path
 */
Path concatenatePaths(const Path& first_path, const Path& second_path);
}
