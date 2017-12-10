//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "graphs/kmer_index.h"

#include <iostream>
#include <list>
#include <string>
#include <unordered_set>
#include <vector>

#include <boost/algorithm/string/join.hpp>

using std::list;
using std::string;
using std::unordered_set;
using std::vector;

struct KmerIndex::KmerIndexImpl {
  explicit KmerIndexImpl(const StringToPathsMap& kmer_to_paths_map);
  explicit KmerIndexImpl(std::shared_ptr<Graph> graph_ptr, int32_t kmer_len);
  void addKmerPathsStartingAtNode(std::shared_ptr<Graph> graph_ptr,
                                  int32_t node_id);
  void addKmerPaths(const std::list<GraphPath>& kmer_paths);
  int32_t _kmer_len;
  StringToPathsMap _kmer_to_paths_map;
};

KmerIndex::KmerIndexImpl::KmerIndexImpl(
    const StringToPathsMap& kmer_to_paths_map)
    : _kmer_to_paths_map(kmer_to_paths_map) {
  _kmer_len = 0;
  for (const auto& kv : _kmer_to_paths_map) {
    const string& kmer = kv.first;
    _kmer_len = kmer.length();
    break;
  }
}

KmerIndex::KmerIndexImpl::KmerIndexImpl(std::shared_ptr<Graph> graph_ptr,
                                        int32_t kmer_len)
    : _kmer_len(kmer_len) {
  for (int32_t node_id = 0; node_id != graph_ptr->NumNodes(); ++node_id) {
    addKmerPathsStartingAtNode(graph_ptr, node_id);
  }
}

void KmerIndex::KmerIndexImpl::addKmerPathsStartingAtNode(
    std::shared_ptr<Graph> graph_ptr, int32_t node_id) {
  const string& node_seq = graph_ptr->NodeSeq(node_id);
  list<int32_t> node_list;
  node_list.push_back(node_id);
  for (size_t pos = 0; pos != node_seq.length(); ++pos) {
    GraphPath path(graph_ptr, pos, node_list, pos);
    addKmerPaths(path.extendBy(0, _kmer_len - 1));
  }
}

void KmerIndex::KmerIndexImpl::addKmerPaths(const list<GraphPath>& kmer_paths) {
  for (const GraphPath& kmer_path : kmer_paths) {
    _kmer_to_paths_map[kmer_path.seq()].push_back(kmer_path);
  }
}

KmerIndex::KmerIndex(const StringToPathsMap& kmer_to_paths_map)
    : _impl(new KmerIndexImpl(kmer_to_paths_map)) {}

KmerIndex::KmerIndex(std::shared_ptr<Graph> graph_ptr, int32_t kmer_len)
    : _impl(new KmerIndexImpl(graph_ptr, kmer_len)) {}

KmerIndex::KmerIndex(const KmerIndex& other)
    : _impl(new KmerIndexImpl(*other._impl)) {}

KmerIndex::KmerIndex(KmerIndex&& other) noexcept
    : _impl(std::move(other._impl)) {}

KmerIndex& KmerIndex::operator=(const KmerIndex& other) {
  if (this != &other) {
    _impl.reset(new KmerIndexImpl(*other._impl));
  }
  return *this;
}

KmerIndex& KmerIndex::operator=(KmerIndex&& other) noexcept {
  _impl = std::move(other._impl);
  return *this;
}

KmerIndex::~KmerIndex() = default;

bool KmerIndex::operator==(const KmerIndex& other) const {
  return (_impl->_kmer_to_paths_map == other._impl->_kmer_to_paths_map &&
          _impl->_kmer_len == other._impl->_kmer_len);
}

static string encodePaths(const list<GraphPath>& paths) {
  list<string> path_encodings;
  for (const auto& path : paths) {
    path_encodings.push_back(path.encode());
  }
  return boost::algorithm::join(path_encodings, ",");
}

string KmerIndex::encode() const {
  list<string> kv_encodings;
  for (const auto& kv : _impl->_kmer_to_paths_map) {
    const string encoding_of_paths = encodePaths(kv.second);
    const string kv_encoding = "{" + kv.first + "->" + encoding_of_paths + "}";
    kv_encodings.push_back(kv_encoding);
  }
  return boost::algorithm::join(kv_encodings, ",");
}

bool KmerIndex::contains(const std::string& kmer) const {
  return _impl->_kmer_to_paths_map.find(kmer) !=
         _impl->_kmer_to_paths_map.end();
}

size_t KmerIndex::numPaths(const std::string& kmer) const {
  if (!contains(kmer)) {
    return 0;
  }
  return _impl->_kmer_to_paths_map.at(kmer).size();
}

const list<GraphPath>& KmerIndex::getPaths(const std::string& kmer) const {
  return _impl->_kmer_to_paths_map.at(kmer);
}

unordered_set<string> KmerIndex::getKmersWithNonzeroCount() const {
  unordered_set<string> kmers_with_nonzero_count;
  for (const auto& kv : _impl->_kmer_to_paths_map) {
    kmers_with_nonzero_count.insert(kv.first);
  }
  return kmers_with_nonzero_count;
}

std::ostream& operator<<(std::ostream& os, const KmerIndex& kmer_index) {
  os << kmer_index.encode();
  return os;
}
