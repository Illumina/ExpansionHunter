//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>
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

#include "graphalign/KmerIndexOperations.hh"

#include <list>

#include <boost/algorithm/string.hpp>

#include "graphutils/SequenceOperations.hh"

using std::list;
using std::string;

namespace graphtools
{
list<string> extractKmersFromAllPositions(const string& sequence, int32_t kmer_len)
{
    list<string> kmers;
    for (size_t pos = 0; pos + kmer_len <= sequence.length(); ++pos)
    {
        string kmer = sequence.substr(pos, static_cast<std::size_t>(kmer_len));
        boost::to_upper(kmer);
        kmers.push_back(kmer);
    }
    return kmers;
}

int32_t countKmerMatches(const KmerIndex& kmer_index, const std::string& seq)
{
    const list<string> kmers = extractKmersFromAllPositions(seq, kmer_index.kmerLength());
    int32_t num_kmer_matches = 0;

    for (const string& kmer : kmers)
    {
        if (kmer_index.numPaths(kmer) != 0)
        {
            ++num_kmer_matches;
        }
    }
    return num_kmer_matches;
}

bool checkIfForwardOriented(const KmerIndex& kmer_index, const std::string& sequence)
{
    const int32_t num_forward_matches = countKmerMatches(kmer_index, sequence);
    const int32_t num_revcomp_matches = countKmerMatches(kmer_index, reverseComplement(sequence));
    return num_forward_matches >= num_revcomp_matches;
}

/**
 * Find minimum kmer length that covers each node with a unique kmer
 * @param graph  a graph
 * @param min_unique_kmers_per_edge  min number of unique kmers to cover each edge
 * @param min_unique_kmers_per_node  min number of unique kmers to cover each node
 * @return
 */
int findMinCoveringKmerLength(Graph const* graph, size_t min_unique_kmers_per_edge, size_t min_unique_kmers_per_node)
{
    for (int32_t k = 10; k < 64; ++k)
    {
        KmerIndex index(*graph, k);

        bool any_below = false;
        for (NodeId node_id = 0; node_id != graph->numNodes(); ++node_id)
        {
            if (index.numUniqueKmersOverlappingNode(node_id) < min_unique_kmers_per_node)
            {
                any_below = true;
                break;
            }
            // this will enumerate all edges
            for (const auto succ : graph->successors(node_id))
            {
                if (index.numUniqueKmersOverlappingEdge(node_id, succ) < min_unique_kmers_per_edge)
                {
                    any_below = true;
                    break;
                }
            }
            if (any_below)
            {
                break;
            }
        }
        if (any_below)
        {
            continue;
        }

        return k;
    }
    return -1;
}
}
