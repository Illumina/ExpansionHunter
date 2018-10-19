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

#include "graphalign/GaplessAligner.hh"

#include <iostream>
#include <stdexcept>
#include <vector>

#include "graphcore/PathOperations.hh"
#include "graphutils/BaseMatching.hh"
#include "graphutils/SequenceOperations.hh"

using std::list;
using std::string;
using std::vector;

namespace graphtools
{
list<GraphAlignment> GaplessAligner::align(const string& query) const
{
    const list<string> kmers = extractKmersFromAllPositions(query, kmer_len_);

    int32_t pos = 0;
    for (const string& kmer : kmers)
    {
        // Initiate alignment from a unique kmer.
        if (kmer_index_.numPaths(kmer) == 1)
        {
            const Path kmer_path = kmer_index_.getPaths(kmer).front();
            return getBestAlignmentToShortPath(kmer_path, pos, query);
        }
        ++pos;
    }
    return {};
}

list<GraphAlignment> getBestAlignmentToShortPath(const Path& path, int32_t start_pos, const string& query)
{
    const int32_t start_extension = start_pos;
    const auto end_extension = static_cast<const int32_t>(query.length() - start_pos - path.length());
    const list<Path> full_paths = extendPath(path, start_extension, end_extension);

    list<GraphAlignment> best_alignments;
    int32_t max_matches = -1;

    for (const Path& full_path : full_paths)
    {
        GraphAlignment alignment = alignWithoutGaps(full_path, query);
        if (static_cast<int32_t>(alignment.numMatches()) > max_matches)
        {
            max_matches = alignment.numMatches();
            best_alignments = { alignment };
        }
        else if (static_cast<int32_t>(alignment.numMatches()) == max_matches)
        {
            best_alignments.push_back(alignment);
        }
    }

    return best_alignments;
}

GraphAlignment alignWithoutGaps(const Path& path, const string& query)
{
    vector<string> query_pieces = splitSequenceByPath(path, query);
    vector<Alignment> alignments;

    const Graph* graph_ptr = path.graphRawPtr();
    size_t index = 0;
    for (auto node_id : path.nodeIds())
    {
        const string& node_seq = graph_ptr->nodeSeq(node_id);
        const string query_piece = query_pieces[index];
        const int32_t ref_start = index == 0 ? path.startPosition() : 0;
        alignments.push_back(alignWithoutGaps(ref_start, node_seq, query_piece));
        ++index;
    }

    return GraphAlignment(path, alignments);
}

Alignment alignWithoutGaps(int32_t ref_start, const string& reference, const string& query)
{
    if (reference.length() < ref_start + query.length())
    {
        throw std::logic_error(
            "Gapless alignment requires that sequences " + query + " and " + reference + " have same length.");
    }

    if (query.empty() || reference.empty())
    {
        throw std::logic_error("Cannot align empty sequences");
    }

    list<Operation> operations;
    size_t previous_run_end = 0;
    size_t run_length = 0;
    char run_operation = '\0';
    for (size_t index = 0; index != query.length(); ++index)
    {
        char cur_operation = 'X';
        if (checkIfReferenceBaseMatchesQueryBase(reference[ref_start + index], query[index]))
        {
            cur_operation = 'M';
        }

        if (cur_operation == run_operation)
        {
            ++run_length;
        }
        else
        {
            if (run_operation != '\0')
            {
                OperationType operation_type = decodeOperationType(run_operation);
                operations.emplace_back(operation_type, run_length);
            }
            previous_run_end += run_length;
            run_length = 1;
            run_operation = cur_operation;
        }
    }

    OperationType operation_type = decodeOperationType(run_operation);
    operations.emplace_back(operation_type, run_length);

    return Alignment(ref_start, operations);
}
}
