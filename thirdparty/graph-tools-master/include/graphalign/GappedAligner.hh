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

#include <cstdint>
#include <list>
#include <string>
#include <utility>

#include "graphalign/GaplessAligner.hh"
#include "graphalign/GraphAligner.hh"
#include "graphalign/GraphAlignment.hh"
#include "graphalign/KmerIndex.hh"
#include "graphalign/LinearAlignment.hh"
#include "graphalign/LinearAlignmentParameters.hh"
#include "graphalign/PinnedDagAligner.hh"
#include "graphalign/PinnedPathAligner.hh"
#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

namespace graphtools
{

using PathAndAlignment = std::pair<Path, Alignment>;

/**
 * General graph aligner supporting linear gaps.
 */
class GappedGraphAligner : public GraphAligner
{
public:
    /**
     * Initializes the aligner
     *
     * @param graph: A graph possibly containing loops (but no cycles)
     * @param kmer_len: Kmer length for kmer index
     * @param padding_len: Elongate paths by this much during path kmer extension step to allow for gaps
     * @param seed_affix_trim_len: Trim length for the prefix and suffix (=affix) of the path
     * @param match_score: Score for matching bases
     * @param mismatch_score: Score for mismatching bases
     * @param gap_open_score: Score for opeaning a gap (linear)
     * @param gap_extend_score: Score for extending an open gap (linear)
     */
    GappedGraphAligner(
        const Graph* graph_ptr, size_t kmer_len, size_t padding_len, size_t seed_affix_trim_len,
        const std::string& alignerName, LinearAlignmentParameters alignerParameters = LinearAlignmentParameters())
        : kmer_len_(kmer_len)
        , padding_len_(padding_len)
        , seed_affix_trim_len_(seed_affix_trim_len)
        , kmer_index_(*graph_ptr, kmer_len)
        , aligner_(alignerName, alignerParameters)
    {
    }

    /**
     * Aligns a read to the graph
     *
     * @param query: Query sequence
     * @return List of top-scoring graph alignments
     */
    std::list<GraphAlignment> align(const std::string& query) const override;

    /**
     * Extends a path matching a kmer in the query sequence to full-length alignments
     *
     * @param kmer_path: Kmer match path
     * @param query: Query sequence
     * @param kmer_start_on_query: Position of the left-most base of the kmer on the query sequence
     * @return List of top-scoring graph alignments going through the kmer match path
     */
    std::list<GraphAlignment>
    extendKmerMatchToFullAlignments(Path kmer_path, const std::string& query, size_t kmer_start_on_query) const;

    /**
     * Aligns query suffix to all suffix-extensions of a given path
     *
     * @param query_piece: Query suffix to align
     * @param seed_path: Path from whose suffix the alignments should start
     * @param extension_len: Length of suffix-extensions
     * @return List of top-scoring alignments and their paths; each path is extended to contain the seed path
     */
    std::list<PathAndAlignment>
    extendAlignmentPrefix(const Path& seed_path, const std::string& query_piece, size_t extension_len) const;

    /**
     * Aligns query prefix to all prefix-extensions of a given path
     *
     * @param query_piece: Query prefix to align
     * @param seed_path: Path at whose prefix the alignments should end
     * @param extension_len: Length of prefix-extensions
     * @return List of top-scoring alignments and their paths; each path is extended to contain the seed path
     */
    std::list<PathAndAlignment>
    extendAlignmentSuffix(const Path& seed_path, const std::string& query_piece, size_t extension_len) const;

private:
    const size_t kmer_len_;
    const size_t padding_len_;
    const int32_t seed_affix_trim_len_;
    const KmerIndex kmer_index_;

    class AlignerSelector
    {
        std::shared_ptr<PinnedPathAligner> ptrPathAligner_;
        mutable std::shared_ptr<PinnedDagAligner> ptrDagAligner_;
    public:
        AlignerSelector(const std::string& alignerName, const LinearAlignmentParameters& alignerParameters)
        {
            if ("path-aligner" == alignerName)
            {
                ptrPathAligner_ = std::make_shared<PinnedPathAligner>(
                    alignerParameters.matchScore, alignerParameters.mismatchScore, alignerParameters.gapOpenScore);
            }
            else if ("dag-aligner" == alignerName)
            {
                ptrDagAligner_ = std::make_shared<PinnedDagAligner>(
                    alignerParameters.matchScore, alignerParameters.mismatchScore, alignerParameters.gapOpenScore,
                    alignerParameters.gapExtendScore);
            }
            else
            {
                throw std::invalid_argument("Aligner " + alignerName + " is not available");
            }
        }

        std::list<PathAndAlignment>
        suffixAlign(const Path& seed_path, const std::string& query_piece, size_t extension_len, int& score) const
        {
            return ptrPathAligner_ ? ptrPathAligner_->suffixAlign(seed_path, query_piece, extension_len, score)
                                   : ptrDagAligner_->suffixAlign(seed_path, query_piece, extension_len, score);
        }

        std::list<PathAndAlignment>
        prefixAlign(const Path& seed_path, const std::string& query_piece, size_t extension_len, int& score) const
        {
            return ptrPathAligner_ ? ptrPathAligner_->prefixAlign(seed_path, query_piece, extension_len, score)
                                   : ptrDagAligner_->prefixAlign(seed_path, query_piece, extension_len, score);
        }


    } aligner_;
};
}
