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

#include <initializer_list>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <regex>
#include <string>
#include <utility>
#include <vector>

#include "graphs/graph.h"

class Operation {
 public:
  enum class Type {
    kMatch,
    kMismatch,
    kInsertionToRef,
    kDeletionFromRef,
    kSoftClipping,
    kMissingBases
  };

  Operation(char type_encoding, int length, std::string query,
            std::string reference) {
    DecodeOperation(type_encoding, length, std::move(query),
                    std::move(reference));
    Validate();
  }
  Operation(std::string cigar, std::string query, std::string reference);
  Type type() const { return type_; }
  int Length() const { return length_; }
  std::string Query() const { return query_; }
  std::string Reference() const { return reference_; }
  int32_t QuerySpan() const;
  int32_t ReferenceSpan() const;
  bool operator==(const Operation& other) const {
    return type_ == other.type_ && length_ == other.length_ &&
           query_ == other.query_ && reference_ == other.reference_;
  }
  char AsSymbol() const;

 private:
  void DecodeOperation(char op_char, int length, std::string query_seq,
                       std::string reference_seq);
  void Validate() const;
  Type type_;
  int length_;
  std::string query_;
  std::string reference_;
};

std::ostream& operator<<(std::ostream& os, const Operation& operation);

class Mapping {
 public:
  typedef size_t size_type;

  Mapping() = default;
  Mapping(int32_t reference_start, std::vector<Operation> operations)
      : reference_start_(reference_start), operations_(std::move(operations)) {
    UpdateCounts();
  }
  Mapping(int32_t reference_start, const std::string& encoding,
          const std::string& query, const std::string& reference);
  size_type NumOperations() const { return operations_.size(); }
  std::string Query() const;
  std::string Reference() const;
  int32_t QuerySpan() const;
  int32_t ReferenceSpan() const;
  int32_t ReferenceStart() const { return reference_start_; }
  void SetReferenceStart(int32_t reference_start) {
    reference_start_ = reference_start;
  }
  const Operation& operator[](size_type index) const {
    return operations_[index];
  }
  bool operator==(const Mapping& other) const {
    return operations_ == other.operations_ &&
           reference_start_ == other.reference_start_;
  }
  size_t NumMatched() const { return matched_; }
  size_t NumMismatched() const { return mismatched_; }
  size_t NumClipped() const { return clipped_; }
  size_t NumInserted() const { return inserted_; }
  size_t NumDeleted() const { return deleted_; }

 protected:
  void DecodeOperations(int32_t reference_start, const std::string& encoding,
                        const std::string& query, const std::string& reference);

  void UpdateCounts();

 private:
  size_t matched_ = 0;
  size_t mismatched_ = 0;
  size_t clipped_ = 0;
  size_t inserted_ = 0;
  size_t deleted_ = 0;
  size_t missing_ = 0;
  int32_t reference_start_ = 0;
  std::vector<Operation> operations_;
};

std::ostream& operator<<(std::ostream& os, const Mapping& mapping);

struct NodeMapping {
  int32_t node_id;
  Mapping mapping;
  bool operator==(const NodeMapping& other) const {
    return node_id == other.node_id && mapping == other.mapping;
  }
};

class GraphMapping {
 public:
  typedef size_t size_type;
  typedef std::vector<NodeMapping> NodeMappings;
  typedef NodeMappings::const_iterator const_iterator;
  GraphMapping() = default;
  GraphMapping(const std::vector<int32_t>& node_ids,
               const std::vector<Mapping>& mappings) {
    for (size_t index = 0; index != node_ids.size(); ++index) {
      NodeMapping node_mapping;
      node_mapping.node_id = node_ids[index];
      node_mapping.mapping = mappings[index];
      node_mappings_.push_back(node_mapping);
    }
  }
  std::string Query() const;
  std::string Reference() const;
  int32_t QuerySpan() const;
  int32_t ReferenceSpan() const;
  int32_t NumMatches() const;
  bool OverlapsNode(int32_t node_id) const;
  std::list<int32_t> GetIndexesOfNode(int32_t node_id) const;
  const_iterator begin() const { return node_mappings_.begin(); }
  const_iterator end() const { return node_mappings_.end(); }
  const NodeMapping& front() const { return node_mappings_.front(); }
  const NodeMapping& back() const { return node_mappings_.back(); }
  size_type size() const { return node_mappings_.size(); }
  const NodeMapping& operator[](size_t index) const {
    return node_mappings_[index];
  }

  bool operator==(const GraphMapping& other) const {
    return node_mappings_ == other.node_mappings_;
  }
  friend std::ostream& operator<<(std::ostream& os,
                                  const GraphMapping& graph_mapping);

 private:
  NodeMappings node_mappings_;
};

std::ostream& operator<<(std::ostream& os, const GraphMapping& graph_mapping);

typedef std::unique_ptr<GraphMapping> GraphMappingPtr;
