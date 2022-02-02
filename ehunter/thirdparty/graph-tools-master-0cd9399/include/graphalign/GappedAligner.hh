//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Roman Petrovski <RPetrovski@illumina.com>
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

#include <cstdint>
#include <list>
#include <string>
#include <utility>

#include <boost/optional.hpp>

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

enum class AlignerType
{
    PATH_ALIGNER,
    DAG_ALIGNER
};

/// Implements alignment details that are independent of the graph
class AlignerSelector
{
    std::shared_ptr<PinnedPathAligner> ptrPathAligner_;
    mutable std::shared_ptr<PinnedDagAligner> ptrDagAligner_;

public:
    explicit AlignerSelector(
        const AlignerType alignerType, const LinearAlignmentParameters& alignerParameters = LinearAlignmentParameters())
    {
        if (alignerType == AlignerType::PATH_ALIGNER)
        {
            ptrPathAligner_ = std::make_shared<PinnedPathAligner>(
                alignerParameters.matchScore, alignerParameters.mismatchScore, alignerParameters.gapOpenScore);
        }
        else if (alignerType == AlignerType::DAG_ALIGNER)
        {
            ptrDagAligner_ = std::make_shared<PinnedDagAligner>(
                alignerParameters.matchScore, alignerParameters.mismatchScore, alignerParameters.gapOpenScore,
                alignerParameters.gapExtendScore);
        }
        else
        {
            throw std::invalid_argument(
                "AlignerType " + std::to_string(static_cast<int>(alignerType)) + " is not available");
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
};

using PathAndAlignment = std::pair<Path, Alignment>;

/**
 * General graph aligner supporting linear gaps.
 */
class GappedGraphAligner
{
public:
    /**
     * @param graph: A graph possibly containing loops (but no cycles)
     * @param kmer_len: Kmer length for kmer index
     * @param padding_len: Elongate paths by this much during path kmer extension step to allow for gaps
     * @param seed_affix_trim_len: Trim length for the prefix and suffix (=affix) of the path
     */
    GappedGraphAligner(const Graph* graph_ptr, size_t kmer_len, size_t padding_len, size_t seed_affix_trim_len)
        : kmer_len_(kmer_len)
        , padding_len_(padding_len)
        , seed_affix_trim_len_(seed_affix_trim_len)
        , kmer_index_(*graph_ptr, kmer_len)
    {
    }

    /**
     * Aligns a read to the graph
     *
     * @param query: Query sequence
     * @return List of top-scoring graph alignments
     */
    std::list<GraphAlignment> align(const std::string& query, AlignerSelector& alignerSelector) const;

    /**
     * Extends a seed path corresponding to a perfect match to the query sequence to full-length alignments
     *
     * @param seed_path: Seed path
     * @param query: Query sequence
     * @param seed_start_on_query: Position of the left-most base of the seed on the query sequence
     * @return List of top-scoring graph alignments going through the seed path
     */
    std::list<GraphAlignment> extendSeedToFullAlignments(
        Path seed_path, const std::string& query, size_t seed_start_on_query, AlignerSelector& alignerSelector) const;

    /**
     * Aligns query suffix to all suffix-extensions of a given path
     *
     * @param query_piece: Query suffix to align
     * @param seed_path: Path from whose suffix the alignments should start
     * @param extension_len: Length of suffix-extensions
     * @return List of top-scoring alignments and their paths; each path is extended to contain the seed path
     */
    std::list<PathAndAlignment> extendAlignmentPrefix(
        const Path& seed_path, const std::string& query_piece, size_t extension_len,
        AlignerSelector& alignerSelector) const;

    /**
     * Aligns query prefix to all prefix-extensions of a given path
     *
     * @param query_piece: Query prefix to align
     * @param seed_path: Path at whose prefix the alignments should end
     * @param extension_len: Length of prefix-extensions
     * @return List of top-scoring alignments and their paths; each path is extended to contain the seed path
     */
    std::list<PathAndAlignment> extendAlignmentSuffix(
        const Path& seed_path, const std::string& query_piece, size_t extension_len,
        AlignerSelector& alignerSelector) const;

private:
    const size_t kmer_len_;
    const size_t padding_len_;
    const int32_t seed_affix_trim_len_;
    const KmerIndex kmer_index_;

    // An alignment seed is a path whose sequence is a perfect match to the query starting from a given position
    struct AlignmentSeed
    {
        AlignmentSeed(Path path, int start_on_query)
            : path(std::move(path))
            , start_on_query(start_on_query)
        {
        }
        Path path;
        int start_on_query = -1;
    };

    // Performs a search for an alignment seed
    boost::optional<AlignmentSeed> searchForAlignmentSeed(const std::string& query) const;
};
}
