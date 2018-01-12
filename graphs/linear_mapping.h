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

#include <string>
#include <vector>

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
  std::string GetCigarString() const;

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
  std::string GetCigarString() const;

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