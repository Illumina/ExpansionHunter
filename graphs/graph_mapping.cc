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

#include <stdexcept>
#include <typeindex>

#include "graphs/graph_mapping.h"

using std::list;
using std::map;
using std::string;
using std::to_string;

Operation::Operation(std::string cigar, string query, string reference) {
  string length_encoding = cigar;
  length_encoding.pop_back();
  char type_encoding = cigar.back();

  DecodeOperation(type_encoding, std::stoi(length_encoding), std::move(query),
                  std::move(reference));
  Validate();
}

void Operation::DecodeOperation(char type_encoding, int length, string query,
                                string reference) {
  query_ = std::move(query);
  reference_ = std::move(reference);

  switch (type_encoding) {
    case 'M':
      type_ = Type::kMatch;
      break;
    case 'N':
      type_ = Type::kMissingBases;
      break;
    case 'X':
      type_ = Type::kMismatch;
      break;
    case 'I':
      type_ = Type::kInsertionToRef;
      break;
    case 'D':
      type_ = Type::kDeletionFromRef;
      break;
    case 'S':
      type_ = Type::kSoftClipping;
      break;
    default:
      throw std::logic_error(to_string(type_encoding) +
                             " is unknown CIGAR operation");
  }
  length_ = length;
}

void Operation::Validate() const {
  const bool full_length_query = query_.length() == (size_t)length_;
  const bool full_length_ref = reference_.length() == (size_t)length_;
  switch (type_) {
    case Type::kMatch:
      if (full_length_query && query_ == reference_) return;
      break;
    case Type::kMismatch:
      if (full_length_query && query_.length() == reference_.length()) {
        bool found_matching_base = false;
        for (size_t index = 0; index != query_.length(); ++index) {
          if (query_[index] == reference_[index]) found_matching_base = true;
        }
        if (!found_matching_base) return;
      }
      break;
    case Type::kMissingBases:
      if (query_.length() == reference_.length() && full_length_query) {
        bool found_non_n_base_in_query_and_ref = false;
        for (size_t index = 0; index != query_.length(); ++index) {
          if (query_[index] != 'N' && reference_[index] != 'N') {
            found_non_n_base_in_query_and_ref = true;
          }
        }
        if (!found_non_n_base_in_query_and_ref) return;
      }
      break;
    case Type::kDeletionFromRef:
      if (query_.empty() && !reference_.empty() && full_length_ref) return;
      break;
    case Type::kInsertionToRef:
      if (!query_.empty() && reference_.empty() && full_length_query) return;
      break;
    case Type::kSoftClipping:
      if (!query_.empty() && reference_.empty() && full_length_query) return;
      break;
  }

  throw std::logic_error(query_ + " and " + reference_ +
                         " are incompatible with operation " + AsSymbol());
}

int32_t Operation::QuerySpan() const {
  switch (type_) {
    case Type::kMatch:
    case Type::kMismatch:
    case Type::kMissingBases:
    case Type::kInsertionToRef:
    case Type::kSoftClipping:
      return length_;
    case Type::kDeletionFromRef:
      return 0;
  }
  return -1;
}

int32_t Operation::ReferenceSpan() const {
  switch (type_) {
    case Type::kMatch:
    case Type::kMismatch:
    case Type::kMissingBases:
    case Type::kDeletionFromRef:
      return length_;
    case Type::kInsertionToRef:
    case Type::kSoftClipping:
      return 0;
  }
  return -1;
}

char Operation::AsSymbol() const {
  static const map<Operation::Type, char> kOperationToChar = {
      {Operation::Type::kMatch, 'M'},
      {Operation::Type::kMismatch, 'X'},
      {Operation::Type::kInsertionToRef, 'I'},
      {Operation::Type::kDeletionFromRef, 'D'},
      {Operation::Type::kSoftClipping, 'S'},
      {Operation::Type::kMissingBases, 'N'}};
  return kOperationToChar.at(type_);
}

std::ostream& operator<<(std::ostream& os, const Operation& operation) {
  os << operation.Length() << operation.AsSymbol() << "("
     << operation.Reference() << "->" << operation.Query() << ")";
  return os;
}

Mapping::Mapping(int32_t reference_start, const std::string& cigar,
                 const std::string& query, const std::string& reference)
    : reference_start_(reference_start) {
  DecodeOperations(reference_start, cigar, query, reference);
  UpdateCounts();
}

void Mapping::UpdateCounts() {
  clipped_ = 0;
  matched_ = 0;
  mismatched_ = 0;
  missing_ = 0;
  inserted_ = 0;
  deleted_ = 0;
  for (size_t i = 0; i < NumOperations(); ++i) {
    auto const& m = operations_[i];
    switch (m.type()) {
      case Operation::Type::kSoftClipping:
        clipped_ += m.Length();
        break;
      case Operation::Type::kMatch:
        matched_ += m.Length();
        break;
      case Operation::Type::kMismatch:
        mismatched_ += m.Length();
        break;
      case Operation::Type::kMissingBases:
        missing_ += m.Length();
        break;
      case Operation::Type::kInsertionToRef:
        inserted_ += m.Length();
        break;
      case Operation::Type::kDeletionFromRef:
        deleted_ += m.Length();
        break;
    }
  }
}

void Mapping::DecodeOperations(int32_t reference_start,
                               const std::string& cigar,
                               const std::string& query,
                               const std::string& reference) {
  int32_t ref_pos = reference_start;
  int32_t query_pos = 0;
  string length_encoding;
  for (char c : cigar) {
    if (isalpha(c) != 0) {
      string query_piece;
      string reference_piece;
      int32_t operation_length = std::stoi(length_encoding);
      switch (c) {
        case 'M':
        case 'X':
        case 'N': {
          query_piece =
              query.substr((size_t)query_pos, (size_t)operation_length);
          reference_piece =
              reference.substr((size_t)ref_pos, (size_t)operation_length);
          query_pos += operation_length;
          ref_pos += operation_length;
          break;
        }
        case 'I':
        case 'S':
          query_piece =
              query.substr((size_t)query_pos, (size_t)operation_length);
          query_pos += operation_length;
          break;
        case 'D':
          reference_piece =
              reference.substr((size_t)ref_pos, (size_t)operation_length);
          ref_pos += operation_length;
          break;
        default:
          throw std::logic_error(to_string(c) + " is unknown CIGAR operation");
      }
      operations_.emplace_back(c, operation_length, query_piece,
                               reference_piece);
      length_encoding.clear();
    } else {
      if (isdigit(c) == 0) {
        throw std::logic_error(cigar + " is malformed CIGAR string");
      }
      length_encoding += c;
    }
  }
}

string Mapping::Query() const {
  string query;
  for (const auto& operation : operations_) {
    if (operation.type() != Operation::Type::kSoftClipping) {
      query += operation.Query();
    }
  }
  return query;
}

string Mapping::Reference() const {
  string reference;
  for (const auto& operation : operations_) {
    reference += operation.Reference();
  }
  return reference;
}

int32_t Mapping::QuerySpan() const {
  int32_t query_span = 0;
  for (const auto& operation : operations_) {
    query_span += operation.QuerySpan();
  }
  return query_span;
}

int32_t Mapping::ReferenceSpan() const {
  int32_t reference_span = 0;
  for (const auto& operation : operations_) {
    reference_span += operation.ReferenceSpan();
  }
  return reference_span;
}

std::ostream& operator<<(std::ostream& os, const Mapping& mapping) {
  os << "Ref start: " << mapping.ReferenceStart() << ", ";
  for (size_t index = 0; index != mapping.NumOperations(); ++index) {
    os << mapping[index];
  }

  return os;
}

string GraphMapping::Query() const {
  string query;
  for (const auto& node_mapping : node_mappings_) {
    query += node_mapping.mapping.Query();
  }
  return query;
}

string GraphMapping::Reference() const {
  string reference;
  for (const auto& node_mapping : node_mappings_) {
    reference += node_mapping.mapping.Reference();
  }
  return reference;
}

int32_t GraphMapping::QuerySpan() const {
  int32_t query_span = 0;
  for (const auto& node_mapping : node_mappings_) {
    query_span += node_mapping.mapping.QuerySpan();
  }
  return query_span;
}

int32_t GraphMapping::ReferenceSpan() const {
  int32_t reference_span = 0;
  for (const auto& node_mapping : node_mappings_) {
    reference_span += node_mapping.mapping.ReferenceSpan();
  }
  return reference_span;
}

int32_t GraphMapping::NumMatches() const {
  int32_t num_matches = 0;
  for (const auto& node_mapping : node_mappings_) {
    num_matches += node_mapping.mapping.NumMatched();
  }
  return num_matches;
}

bool GraphMapping::OverlapsNode(int32_t node_id) const {
  for (const auto& node_mapping : node_mappings_) {
    if (node_mapping.node_id == node_id) {
      return true;
    }
  }
  return false;
}

list<int32_t> GraphMapping::GetIndexesOfNode(int32_t node_id) const {
  list<int32_t> indexes;
  for (int32_t node_index = 0; node_index != (int32_t)node_mappings_.size();
       ++node_index) {
    if (node_mappings_[node_index].node_id == node_id) {
      indexes.push_back(node_index);
    }
  }
  return indexes;
}

std::ostream& operator<<(std::ostream& os, const GraphMapping& graph_mapping) {
  for (const NodeMapping& node_mapping : graph_mapping) {
    os << node_mapping.node_id << "[" << node_mapping.mapping << "]";
  }
  return os;
}
