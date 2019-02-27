//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>
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

#include "graphalign/KmerIndex.hh"

#include <iostream>
#include <list>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graphcore/PathOperations.hh"
#include "graphutils/PairHashing.hh"
#include "graphutils/SequenceOperations.hh"

#include <boost/algorithm/string/join.hpp>

using std::list;
using std::string;
using std::unordered_set;
using std::vector;

namespace graphtools
{

struct KmerIndex::Impl
{
    explicit Impl(StringToPathsMap kmer_to_paths_map_);
    explicit Impl(const Graph& graph, size_t kmer_len_);
    void addKmerPathsStartingAtNode(const Graph& graph, NodeId node_id);
    void addKmerPaths(const std::list<Path>& kmer_paths);
    void updateKmerCounts();
    size_t kmer_len;
    StringToPathsMap kmer_to_paths_map;
    std::unordered_map<NodeId, size_t> node_kmer_counts;
    std::unordered_map<NodeIdPair, size_t> edge_kmer_counts;
};

KmerIndex::Impl::Impl(StringToPathsMap kmer_to_paths_map_)
    : kmer_to_paths_map(std::move(kmer_to_paths_map_))
{
    kmer_len = 0;
    for (const auto& kv : kmer_to_paths_map)
    {
        const string& kmer = kv.first;
        kmer_len = kmer.length();
        break;
    }
    updateKmerCounts();
}

KmerIndex::Impl::Impl(const Graph& graph, size_t kmer_len_)
    : kmer_len(kmer_len_)
{
    for (NodeId node_id = 0; node_id != graph.numNodes(); ++node_id)
    {
        addKmerPathsStartingAtNode(graph, node_id);
    }
    updateKmerCounts();
}

void KmerIndex::Impl::addKmerPathsStartingAtNode(const Graph& graph, NodeId node_id)
{
    const string node_seq = graph.nodeSeq(node_id);
    vector<NodeId> node_list;
    node_list.push_back(node_id);
    for (size_t pos = 0; pos != node_seq.length(); ++pos)
    {
        Path path(&graph, static_cast<int32_t>(pos), node_list, static_cast<int32_t>(pos));
        addKmerPaths(extendPath(path, 0, static_cast<int32_t>(kmer_len)));
    }
}

void KmerIndex::Impl::addKmerPaths(const list<Path>& kmer_paths)
{
    for (const Path& kmer_path : kmer_paths)
    {
        vector<string> expanded_sequences;
        if (kmer_path.graphRawPtr()->isSequenceExpansionRequired())
        {
            expandReferenceSequence(kmer_path.seq(), expanded_sequences);
        }
        else
        {
            expanded_sequences = { kmer_path.seq() };
        }
        for (const auto& expanded_kmer_seq : expanded_sequences)
        {
            kmer_to_paths_map[expanded_kmer_seq].push_back(kmer_path);
        }
    }
}

void KmerIndex::Impl::updateKmerCounts()
{
    node_kmer_counts.clear();
    edge_kmer_counts.clear();
    for (const auto& kmer_and_paths : kmer_to_paths_map)
    {
        // kmer is unique
        if (kmer_and_paths.second.size() == 1)
        {
            bool has_previous = false;
            NodeId previous_node = 0;
            for (auto const& path_node_id : kmer_and_paths.second.front().nodeIds())
            {
                node_kmer_counts[path_node_id] += 1;
                if (has_previous)
                {
                    edge_kmer_counts[std::make_pair(previous_node, path_node_id)] += 1;
                }
                has_previous = true;
                previous_node = path_node_id;
            }
        }
    }
}

KmerIndex::KmerIndex(const StringToPathsMap& kmer_to_paths_map)
    : pimpl_(new Impl(kmer_to_paths_map))
{
}

KmerIndex::KmerIndex(const Graph& graph, int32_t kmer_len)
    : pimpl_(new Impl(graph, static_cast<size_t>(kmer_len)))
{
}

KmerIndex::KmerIndex(const KmerIndex& other)
    : pimpl_(new Impl(*other.pimpl_))
{
}

KmerIndex::KmerIndex(KmerIndex&& other) noexcept
    : pimpl_(std::move(other.pimpl_))
{
}

KmerIndex& KmerIndex::operator=(const KmerIndex& other)
{
    if (this != &other)
    {
        pimpl_.reset(new Impl(*other.pimpl_));
    }
    return *this;
}

KmerIndex& KmerIndex::operator=(KmerIndex&& other) noexcept
{
    pimpl_ = std::move(other.pimpl_);
    return *this;
}

KmerIndex::~KmerIndex() = default;

bool KmerIndex::operator==(const KmerIndex& other) const
{
    return (pimpl_->kmer_to_paths_map == other.pimpl_->kmer_to_paths_map && pimpl_->kmer_len == other.pimpl_->kmer_len);
}

static string encodePaths(const list<Path>& paths)
{
    list<string> path_encodings;
    for (const auto& path : paths)
    {
        path_encodings.push_back(path.encode());
    }
    return boost::algorithm::join(path_encodings, ",");
}

size_t KmerIndex::kmerLength() const { return pimpl_->kmer_len; }

string KmerIndex::encode() const
{
    list<string> kv_encodings;
    for (const auto& kv : pimpl_->kmer_to_paths_map)
    {
        const string encoding_of_paths = encodePaths(kv.second);
        const string kv_encoding = "{" + kv.first + "->" + encoding_of_paths + "}";
        kv_encodings.push_back(kv_encoding);
    }
    return boost::algorithm::join(kv_encodings, ",");
}

bool KmerIndex::contains(const std::string& kmer) const
{
    return pimpl_->kmer_to_paths_map.find(kmer) != pimpl_->kmer_to_paths_map.end();
}

size_t KmerIndex::numPaths(const std::string& kmer) const
{
    if (!contains(kmer))
    {
        return 0;
    }
    return pimpl_->kmer_to_paths_map.at(kmer).size();
}

const list<Path>& KmerIndex::getPaths(const std::string& kmer) const { return pimpl_->kmer_to_paths_map.at(kmer); }

unordered_set<string> KmerIndex::kmers() const
{
    unordered_set<string> kmers;
    for (const auto& kv : pimpl_->kmer_to_paths_map)
    {
        kmers.insert(kv.first);
    }
    return kmers;
}

size_t KmerIndex::numUniqueKmersOverlappingNode(NodeId node_id) const
{
    auto node_it = pimpl_->node_kmer_counts.find(node_id);
    if (node_it != pimpl_->node_kmer_counts.end())
    {
        return node_it->second;
    }
    return 0;
}

size_t KmerIndex::numUniqueKmersOverlappingEdge(NodeId from, NodeId to) const
{
    auto edge_it = pimpl_->edge_kmer_counts.find(std::make_pair(from, to));
    if (edge_it != pimpl_->edge_kmer_counts.end())
    {
        return edge_it->second;
    }
    return 0;
}

std::ostream& operator<<(std::ostream& os, const KmerIndex& kmer_index)
{
    os << kmer_index.encode();
    return os;
}
}
