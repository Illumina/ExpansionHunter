//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Peter Krusche <pkrusche@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>
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

#include <string>

#include "graphalign/KmerIndex.hh"

namespace graphtools
{

// Returns true if the sequence is in graph's (forward) orientation
bool checkIfForwardOriented(const KmerIndex& kmer_index, const std::string& sequence);

/**
 * Find minimum kmer length that covers each node with a unique kmer
 * @param graph  a graph
 * @param min_unique_kmers_per_edge  min number of unique kmers to cover each edge
 * @param min_unique_kmers_per_node  min number of unique kmers to cover each node
 * @return
 */
int findMinCoveringKmerLength(Graph const* graph, size_t min_unique_kmers_per_edge, size_t min_unique_kmers_per_node);
}
