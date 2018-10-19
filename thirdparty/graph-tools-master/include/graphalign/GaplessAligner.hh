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
#include <string>

#include "graphalign/GraphAligner.hh"
#include "graphalign/GraphAlignment.hh"
#include "graphalign/KmerIndex.hh"
#include "graphcore/Path.hh"

namespace graphtools
{
/**
 * Gapless graph aligner
 *
 * Enables gapless alignment of any sequence to a graph.
 */
class GaplessAligner : public GraphAligner
{
public:
    GaplessAligner(const Graph* graph_ptr, int32_t kmer_len)
        : kmer_len_(kmer_len)
        , kmer_index_(*graph_ptr, kmer_len)
    {
    }
    std::list<GraphAlignment> align(const std::string& query) const override;

private:
    int32_t kmer_len_;
    KmerIndex kmer_index_;
};

/**
 * Determines orientation of a read relative to the graph
 */
class StrandClassifier
{
public:
    StrandClassifier(const Graph& graph, int32_t kmer_len)
        : kmer_len_(kmer_len)
        , kmer_index_(graph, kmer_len)
    {
    }
    bool isForwardOriented(const std::string& seq) const;

private:
    int32_t countKmerMatches(const std::string& seq) const;
    int32_t kmer_len_;
    KmerIndex kmer_index_;
};

/**
 * Computes a top-scoring gapless alignment of a query sequence to the graph
 * that goes through the path starting at the given position on the sequence
 *
 * @param path: Any path shorter than the query
 * @param start_pos: Position on the query corrsponding to the start of the
 * path
 * @param query: Any sequence
 * @return Best gapless alignment with the above properety
 */
std::list<GraphAlignment> getBestAlignmentToShortPath(const Path& path, int32_t start_pos, const std::string& query);

/**
 * Aligns a query sequence to a path of the same length
 *
 * @param path: Any graph path
 * @param query: Any sequence that has the same length as the path
 * @return Result of the alignment
 */
GraphAlignment alignWithoutGaps(const Path& path, const std::string& query);

/**
 * Aligns query sequence to the reference sequence starting at the given
 * position
 *
 * @param query: Any sequence
 * @param ref_start: Position of the start of the alignment on the reference
 * @param reference: Any sequence
 * @return Result of the alignment
 */
Alignment alignWithoutGaps(int32_t ref_start, const std::string& reference, const std::string& query);

/**
 * Extracts kmers starting at each position
 *
 * @param query: Any sequence
 * @param kmer_len: Kmer length
 * @return List of kmers indexed by start position in the original sequence
 */
std::list<std::string> extractKmersFromAllPositions(const std::string& query, int32_t kmer_len);
}
