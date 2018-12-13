//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>

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
