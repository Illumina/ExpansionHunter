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
