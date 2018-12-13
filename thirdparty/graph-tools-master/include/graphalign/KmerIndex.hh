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

#pragma once

#include <list>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

namespace graphtools
{

typedef std::unordered_map<std::string, std::list<Path>> StringToPathsMap;

// Kmer index holds paths that correspond to each kmer that appears in the graph and supports a few standard operations.
class KmerIndex
{
public:
    explicit KmerIndex(const Graph& graph, int32_t kmer_len = 12);
    explicit KmerIndex(const StringToPathsMap& kmer_to_paths_map);
    KmerIndex(const KmerIndex& other);
    KmerIndex(KmerIndex&& other) noexcept;
    KmerIndex& operator=(const KmerIndex& other);
    KmerIndex& operator=(KmerIndex&& other) noexcept;
    ~KmerIndex();
    bool operator==(const KmerIndex& other) const;
    std::string encode() const;
    const std::list<Path>& getPaths(const std::string& kmer) const;
    bool contains(const std::string& kmer) const;
    size_t numPaths(const std::string& kmer) const;
    std::unordered_set<std::string> kmers() const;
    size_t kmerLength() const;

    size_t numUniqueKmersOverlappingNode(NodeId node_id) const;
    size_t numUniqueKmersOverlappingEdge(NodeId from, NodeId to) const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

std::ostream& operator<<(std::ostream& os, const KmerIndex& kmer_index);
}
