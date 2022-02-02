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
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <sparsepp/spp.h>

#include "graphcore/PathOperations.hh"
#include "graphutils/KmerEncoding.hh"
#include "graphutils/PairHashing.hh"
#include "graphutils/SequenceOperations.hh"

using std::list;
using std::string;
using std::unordered_set;
using std::vector;

namespace graphtools
{

struct MiniPath
{
    bool operator==(const MiniPath& other) const
    {
        return (start_position == other.start_position) and (node_id == other.node_id);
    }

    using pos_t = uint16_t;
    using node_t = uint16_t;
    pos_t start_position;
    node_t node_id;
};

using KmerKey_t = uint32_t;

typedef spp::sparse_hash_map<KmerKey_t, MiniPath> StringToMiniPathsMap;
typedef spp::sparse_hash_map<KmerKey_t, std::vector<Path>> StringToPathsMap;

struct KmerIndex::Impl
{
    explicit Impl(const Graph& graph, size_t kmer_len_);
    void addKmerPathsStartingAtNode(NodeId node_id);
    void addKmerPaths(const std::list<Path>& kmer_paths);
    void updateKmerCounts();
    Path MiniPathToPath(const MiniPath& miniPath) const;

private:
    const Graph& graph_;

public:
    size_t kmer_len;
    TwoBitKmerEncoder kmer_coder;
    StringToMiniPathsMap kmer_to_minipaths_map;
    StringToPathsMap kmer_to_paths_map;
    std::unordered_map<NodeId, size_t> node_kmer_counts;
    std::unordered_map<NodeIdPair, size_t> edge_kmer_counts;
};

KmerIndex::Impl::Impl(const Graph& graph, size_t kmer_len_)
    : graph_(graph)
    , kmer_len(kmer_len_)
    , kmer_coder(kmer_len)
{
    for (NodeId node_id = 0; node_id != graph.numNodes(); ++node_id)
    {
        addKmerPathsStartingAtNode(node_id);
    }
    updateKmerCounts();
}

void KmerIndex::Impl::addKmerPathsStartingAtNode(NodeId node_id)
{
    const string node_seq = graph_.nodeSeq(node_id);
    vector<NodeId> node_list;
    node_list.push_back(node_id);
    for (size_t pos = 0; pos != node_seq.length(); ++pos)
    {
        Path path(&graph_, static_cast<int32_t>(pos), node_list, static_cast<int32_t>(pos));
        addKmerPaths(extendPathEnd(path, static_cast<int32_t>(kmer_len)));
    }
}

Path KmerIndex::Impl::MiniPathToPath(const MiniPath& miniPath) const
{
    return Path(&graph_, miniPath.start_position, { miniPath.node_id }, miniPath.start_position + kmer_len);
}

void KmerIndex::Impl::addKmerPaths(const list<Path>& kmer_paths)
{
    for (const Path& kmer_path : kmer_paths)
    {
        vector<string> expanded_sequences;
        expandReferenceSequence(kmer_path.seq(), expanded_sequences);
        for (const auto& expanded_kmer_seq : expanded_sequences)
        {
            const auto expanded_kmer_key = kmer_coder.encode(expanded_kmer_seq);
            auto pathIter(kmer_to_paths_map.find(expanded_kmer_key));
            if (pathIter != kmer_to_paths_map.end())
            {
                pathIter->second.push_back(kmer_path);
                continue;
            }

            auto miniPathIter(kmer_to_minipaths_map.find(expanded_kmer_key));
            if (miniPathIter != kmer_to_minipaths_map.end())
            {
                kmer_to_paths_map[expanded_kmer_key] = { MiniPathToPath(miniPathIter->second), kmer_path };
                kmer_to_minipaths_map.erase(miniPathIter);
            }
            else
            {
                if ((kmer_path.numNodes() == 1)
                    and ((kmer_path.endPosition() - kmer_path.startPosition()) == static_cast<int>(kmer_len))
                    and (kmer_path.startPosition() >= 0)
                    and (kmer_path.startPosition() <= std::numeric_limits<MiniPath::pos_t>::max())
                    and (kmer_path.nodeIds().front() <= std::numeric_limits<MiniPath::node_t>::max()))
                {
                    MiniPath miniPath{ static_cast<MiniPath::pos_t>(kmer_path.startPosition()),
                                       static_cast<MiniPath::node_t>(kmer_path.nodeIds().front()) };
                    kmer_to_minipaths_map.emplace(expanded_kmer_key, miniPath);
                    continue;
                }
                else
                {
                    kmer_to_paths_map[expanded_kmer_key] = { kmer_path };
                }
            }
        }
    }
}

void KmerIndex::Impl::updateKmerCounts()
{
    node_kmer_counts.clear();
    edge_kmer_counts.clear();
    for (const auto& kmer_and_paths : kmer_to_minipaths_map)
    {
        auto const path_node_id(kmer_and_paths.second.node_id);
        node_kmer_counts[path_node_id] += 1;
    }

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
    return (
        pimpl_->kmer_to_minipaths_map == other.pimpl_->kmer_to_minipaths_map
        and pimpl_->kmer_to_paths_map == other.pimpl_->kmer_to_paths_map && pimpl_->kmer_len == other.pimpl_->kmer_len);
}

static string encodePaths(const std::vector<Path>& paths)
{
    std::vector<string> path_encodings;
    for (const auto& path : paths)
    {
        path_encodings.push_back(path.encode());
    }
    return boost::algorithm::join(path_encodings, ",");
}

size_t KmerIndex::kmerLength() const { return pimpl_->kmer_len; }

string KmerIndex::encode() const
{
    std::vector<string> kv_encodings;
    for (const auto& kv : pimpl_->kmer_to_minipaths_map)
    {
        const string encoding_of_paths = encodePaths({ pimpl_->MiniPathToPath(kv.second) });
        kv_encodings.emplace_back("{" + pimpl_->kmer_coder.decode(kv.first) + "->" + encoding_of_paths + "}");
    }
    for (const auto& kv : pimpl_->kmer_to_paths_map)
    {
        const string encoding_of_paths = encodePaths(kv.second);
        kv_encodings.emplace_back("{" + pimpl_->kmer_coder.decode(kv.first) + "->" + encoding_of_paths + "}");
    }
    return boost::algorithm::join(kv_encodings, ",");
}

bool KmerIndex::contains(const std::string& kmer) const
{
    if ((kmer.size() != pimpl_->kmer_len) or (kmer.find_last_not_of("ACGT") != std::string::npos))
    {
        return false;
    }

    const auto kmer_key = pimpl_->kmer_coder.encode(kmer);
    return (
        (pimpl_->kmer_to_minipaths_map.find(kmer_key) != pimpl_->kmer_to_minipaths_map.end())
        or (pimpl_->kmer_to_paths_map.find(kmer_key) != pimpl_->kmer_to_paths_map.end()));
}

size_t KmerIndex::numPaths(const std::string& kmer) const
{
    if ((kmer.size() != pimpl_->kmer_len) or (kmer.find_last_not_of("ACGT") != std::string::npos))
    {
        return 0;
    }

    const auto kmer_key = pimpl_->kmer_coder.encode(kmer);
    if (pimpl_->kmer_to_minipaths_map.find(kmer_key) != pimpl_->kmer_to_minipaths_map.end())
    {
        return 1;
    }
    else
    {
        const auto pathIter(pimpl_->kmer_to_paths_map.find(kmer_key));
        if (pathIter == pimpl_->kmer_to_paths_map.end())
        {
            return 0;
        }
        else
        {
            return pathIter->second.size();
        }
    }
}

std::vector<Path> KmerIndex::getPaths(const std::string& kmer) const
{
    const auto kmer_key = pimpl_->kmer_coder.encode(kmer);
    const auto miniPathIter(pimpl_->kmer_to_minipaths_map.find(kmer_key));
    if (miniPathIter != pimpl_->kmer_to_minipaths_map.end())
    {
        return { pimpl_->MiniPathToPath(miniPathIter->second) };
    }
    else
    {
        return pimpl_->kmer_to_paths_map.at(kmer_key);
    }
}

unordered_set<string> KmerIndex::kmers() const
{
    unordered_set<string> kmers;
    for (const auto& kv : pimpl_->kmer_to_minipaths_map)
    {
        kmers.insert(pimpl_->kmer_coder.decode(kv.first));
    }
    for (const auto& kv : pimpl_->kmer_to_paths_map)
    {
        kmers.insert(pimpl_->kmer_coder.decode(kv.first));
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
