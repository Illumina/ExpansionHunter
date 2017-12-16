//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#pragma once

#include <list>
#include <memory>
#include <string>

#include "graphs/graph_mapping.h"
#include "graphs/kmer_index.h"
#include "graphs/path.h"

/**
 * @brief Gapless graph aligner
 *
 * Enables gapless alignment of any sequence to a graph.
 */
class GaplessAligner {
 public:
  GaplessAligner(std::shared_ptr<Graph> graph_ptr, int32_t kmer_len)
      : kmer_len_(kmer_len), kmer_index_(graph_ptr, kmer_len) {}
  GraphMapping GetBestAlignment(const std::string& sequence) const;

 private:
  int32_t kmer_len_;
  KmerIndex kmer_index_;
};

/**
 * @brief Determines orientation of a read relative to the graph
 *
 */
class StrandClassifier {
 public:
  StrandClassifier(std::shared_ptr<Graph> graph_ptr, int32_t kmer_len)
      : kmer_len_(kmer_len), kmer_index_(graph_ptr, kmer_len) {}
  bool IsForwardOriented(const std::string& seq) const;

 private:
  int32_t CountKmerMatches(const std::string& seq) const;
  int32_t kmer_len_;
  KmerIndex kmer_index_;
};

/**
 * @brief Computes a top-scoring gapless alignment of a sequence to the graph
 * that goes through the path starting at the given position on the sequence
 *
 * @param path: Any path shorter than the sequence
 * @param start_pos: Position on the sequence corrsponding to the start of the
 * path
 * @param sequence: Any sequence
 * @return Best gapless alignment with the above properety
 */
GraphMapping GetBestAlignmentToShortPath(const GraphPath& path,
                                         int32_t start_pos,
                                         const std::string& sequence);

/**
 * @brief Aligns a sequence to a path of the same length
 *
 * @param path: Any graph path
 * @param sequence: Any sequence that has the same length as the path
 * @return Result of the alignment
 */
GraphMapping AlignWithoutGaps(const GraphPath& path,
                              const std::string& sequence);

/**
 * @brief Aligns query sequence to the reference sequence starting at the given
 * position
 *
 * @param query: Any sequence
 * @param ref_start: Position of the start of the alignment on the reference
 * @param reference: Any sequence
 * @return Result of the alignment
 */
Mapping AlignWithoutGaps(const std::string& query, int32_t ref_start,
                         const std::string& reference);

/**
 * @brief Extracts kmers starting at each position
 *
 * @param sequence: Any sequence
 * @param kmer_len: Kmer length
 * @return List of kmers indexed by start position in the original sequence
 */
std::list<std::string> ExtractKmersFromAllPositions(const std::string& sequence,
                                                    int32_t kmer_len);
