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
