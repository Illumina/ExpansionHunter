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

#pragma once

#include <unordered_map>
#include <unordered_set>
#include <memory>

#include "graphs/path.h"

typedef std::unordered_map<std::string, std::list<GraphPath>> StringToPathsMap;

// Kmer index holds paths that correspond to each kmer that appears in the graph
// and supports a few standard operations.
class KmerIndex {
 public:
  explicit KmerIndex(std::shared_ptr<Graph> wgraph, int32_t kmer_len = 12);
  explicit KmerIndex(const StringToPathsMap& kmer_to_paths_map);
  explicit KmerIndex(const KmerIndex& other);
  explicit KmerIndex(KmerIndex&& other) noexcept;
  KmerIndex& operator=(const KmerIndex& other);
  KmerIndex& operator=(KmerIndex&& other) noexcept;
  ~KmerIndex();
  bool operator==(const KmerIndex& other) const;
  std::string encode() const;
  const std::list<GraphPath>& getPaths(const std::string& kmer) const;
  bool contains(const std::string& kmer) const;
  size_t numPaths(const std::string& kmer) const;
  std::unordered_set<std::string> getKmersWithNonzeroCount() const;

 private:
  struct KmerIndexImpl;
  std::unique_ptr<KmerIndexImpl> _impl;
};

std::ostream& operator<<(std::ostream& os, const KmerIndex& kmer_index);
