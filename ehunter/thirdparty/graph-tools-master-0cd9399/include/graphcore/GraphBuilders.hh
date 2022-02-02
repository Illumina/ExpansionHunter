//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>

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

#include <string>

#include "graphcore/Graph.hh"

namespace graphtools
{

/**
 * Construct a graph representing deletion of a sequence from a reference
 *
 * The node ids are assigned in order specified by the function paramters.
 *
 * @param left_flank Sequence of the left flank
 * @param deletion Sequence deleted from the reference
 * @param right_flank Sequence of the right flank
 * @return Deletion graph
 */
Graph makeDeletionGraph(const std::string& left_flank, const std::string& deletion, const std::string& right_flank);

/**
 * Construct a graph representing replacement of a piece of a reference by another sequence
 *
 * The node ids are assigned in order specified by the function paramters.
 *
 * @param left_flank Sequence of the left flank
 * @param deletion Sequence deleted from the reference
 * @param insertion Sequence inserted into the reference
 * @param right_flank Sequence of the right flank
 * @return Swap graph
 */
Graph makeSwapGraph(
    const std::string& left_flank, const std::string& deletion, const std::string& insertion,
    const std::string& right_flank);

/**
 * Construct a graph representing two sequence swaps separated by another sequence
 *
 * The node ids are assigned in order specified by the function paramters.
 *
 * @param left_flank Sequence of the left flank
 * @param deletion1 Deleted sequence of the first swap
 * @param insertion1 Inserted sequence of the first swap
 * @param middle Sequence separating the swaps
 * @param deletion2 Deleted sequence of the second swap
 * @param insertion2 Inserted sequence of the second swap
 * @param right_flank Sequence of the right flank
 * @return Double-swap graph
 */
Graph makeDoubleSwapGraph(
    const std::string& left_flank, const std::string& deletion1, const std::string& insertion1,
    const std::string& middle, const std::string& deletion2, const std::string& insertion2,
    const std::string& right_flank);

/**
 * Construct a graph representing an STR repeat with the loop separated into multiple nodes to keep the graph acyclic
 *
 * The first and the last nodes correspond to the left and the right flanks respectively. The internal nodes correspond
 * to the repeat unit. The number of repeat unit nodes is given by ceiling(read length/unit length).
 *
 * @param read_len Length of the sequenced reads
 * @param left_flank Sequence of the left flank
 * @param repeat_unit Repeat unit
 * @param right_flank Sequence of the right flank
 * @return Loopless STR graph
 */
Graph makeLooplessStrGraph(
    size_t read_len, const std::string& left_flank, const std::string& repeat_unit, const std::string& right_flank);

/**
 * Construct a graph representing an STR repeat
 *
 * The graph consists of the repeat flanks separated by the loop corresponding to multiple repetitions of the repeat
 * unit. The node ids are assigned in order specified by the function paramters.
 *
 * @param left_flank Sequence of the left flank
 * @param repeat_unit Repeat unit
 * @param right_flank Sequence of the right flank
 * @return STR graph
 */
Graph makeStrGraph(const std::string& left_flank, const std::string& repeat_unit, const std::string& right_flank);
}
