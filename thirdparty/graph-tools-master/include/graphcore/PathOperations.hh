// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.

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
