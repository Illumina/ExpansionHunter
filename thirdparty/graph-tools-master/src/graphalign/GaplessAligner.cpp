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
    const auto end_extension = static_cast<int32_t>(query.length() - start_pos - path.length());
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
