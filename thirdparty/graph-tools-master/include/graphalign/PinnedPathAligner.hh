//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
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

#include "graphalign/GraphAlignment.hh"
#include "graphalign/LinearAlignmentOperations.hh"
#include "graphalign/PinnedAligner.hh"
#include "graphcore/PathOperations.hh"

namespace graphtools
{

using PathAndAlignment = std::pair<Path, Alignment>;

class PinnedPathAligner
{
    const int32_t matchScore_;
    const int32_t mismatchScore_;
    const int32_t gapOpenScore_;

    mutable PinnedAligner pinnedAligner_;

public:
    PinnedPathAligner(int32_t matchScore = 5, int32_t mismatchScore = -4, int32_t gapOpenScore = -8)
        : matchScore_(matchScore)
        , mismatchScore_(mismatchScore)
        , gapOpenScore_(gapOpenScore)
        , pinnedAligner_(matchScore_, mismatchScore_, gapOpenScore_)
    {
    }
    std::list<PathAndAlignment>
    suffixAlign(const Path& seed_path, const std::string& query_piece, size_t extension_len, int& score) const;
    std::list<PathAndAlignment>
    prefixAlign(const Path& seed_path, const std::string& query_piece, size_t extension_len, int& score) const;

private:
    int32_t scoreAlignment(const Alignment& alignment) const
    {
        return graphtools::scoreAlignment(alignment, matchScore_, mismatchScore_, gapOpenScore_);
    }
};

inline std::list<PathAndAlignment> PinnedPathAligner::suffixAlign(
    const Path& seed_path, const std::string& query_piece, size_t extension_len, int& top_alignment_score) const
{
    std::list<PathAndAlignment> top_paths_and_alignments;
    top_alignment_score = INT32_MIN;

    const std::list<Path> path_extensions = extendPathStart(seed_path, extension_len);
    for (const auto& path : path_extensions)
    {
        Alignment alignment = pinnedAligner_.suffixAlign(path.seq(), query_piece);
        const int32_t alignment_score = scoreAlignment(alignment);

        if (top_alignment_score < alignment_score)
        {
            top_paths_and_alignments.clear();
            top_alignment_score = alignment_score;
        }

        if (top_alignment_score == alignment_score)
        {
            top_paths_and_alignments.push_back(std::make_pair(path, alignment));
        }
    }

    return top_paths_and_alignments;
}

inline std::list<PathAndAlignment> PinnedPathAligner::prefixAlign(
    const Path& seed_path, const std::string& query_piece, size_t extension_len, int& top_alignment_score) const
{
    std::list<PathAndAlignment> top_paths_and_alignments;
    top_alignment_score = INT32_MIN;

    const std::list<Path> path_extensions = extendPathEnd(seed_path, extension_len);
    for (const auto& path : path_extensions)
    {
        Alignment alignment = pinnedAligner_.prefixAlign(path.seq(), query_piece);
        const int32_t alignment_score = scoreAlignment(alignment);

        if (top_alignment_score < alignment_score)
        {
            top_paths_and_alignments.clear();
            top_alignment_score = alignment_score;
        }

        if (top_alignment_score == alignment_score)
        {
            top_paths_and_alignments.push_back(std::make_pair(path, alignment));
        }
    }

    return top_paths_and_alignments;
}
}
