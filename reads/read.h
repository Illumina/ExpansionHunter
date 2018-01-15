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

#include <cstdint>
#include <memory>
#include <sstream>
#include <string>

#include "classification/mapping_classifier.h"
#include "graphs/graph_mapping.h"

namespace reads {

struct CoreInfo {
  std::string fragment_id;
  std::string bases;
  std::string quals;
  bool operator==(const CoreInfo& other) const {
    return (fragment_id == other.fragment_id && bases == other.bases &&
            quals == other.quals);
  }
};

struct SamInfo {
  int32_t chrom_id = -1;
  int32_t pos = -1;
  int32_t mapq = -1;
  int32_t mate_chrom_id = -1;
  int32_t mate_pos = -1;
  bool is_mapped = false;
  bool is_first_mate = false;
  bool is_mate_mapped = false;
  bool operator==(const SamInfo& other) const {
    return (chrom_id == other.chrom_id && pos == other.pos &&
            mapq == other.mapq && mate_chrom_id == other.mate_chrom_id &&
            mate_pos == other.mate_pos && is_mapped == other.is_mapped &&
            is_first_mate == other.is_first_mate &&
            is_mate_mapped == other.is_mate_mapped);
  }
};

struct GraphInfo {
  GraphInfo() = default;
  GraphInfo(const GraphInfo& other)
      : canonical_mapping_type(other.canonical_mapping_type),
        num_str_units_spanned(other.num_str_units_spanned) {
    if (other.canonical_mapping_ptr) {
      canonical_mapping_ptr.reset(
          new GraphMapping(*other.canonical_mapping_ptr));
    }
  }
  GraphInfo& operator=(const GraphInfo& other) {
    if (other.canonical_mapping_ptr) {
      canonical_mapping_ptr.reset(
          new GraphMapping(*other.canonical_mapping_ptr));
    }
    canonical_mapping_type = other.canonical_mapping_type;
    num_str_units_spanned = other.num_str_units_spanned;
    return *this;
  }
  std::unique_ptr<GraphMapping> canonical_mapping_ptr;
  MappingType canonical_mapping_type = MappingType::kUnmapped;
  int32_t num_str_units_spanned = 0;
  bool operator==(const GraphInfo& other) const {
    if (canonical_mapping_type != other.canonical_mapping_type ||
        num_str_units_spanned != other.num_str_units_spanned) {
      return false;
    }

    if (!canonical_mapping_ptr && !other.canonical_mapping_ptr) {
      return true;
    }

    if (!canonical_mapping_ptr || !other.canonical_mapping_ptr) {
      return false;
    }

    return (*canonical_mapping_ptr == *other.canonical_mapping_ptr);
  }
};

class Read {
 public:
  Read(const std::string& fragment_id, const std::string& bases,
       const std::string& quals) {
    SetCoreInfo(fragment_id, bases, quals);
  }

  void SetCoreInfo(const std::string& fragment_id, const std::string& bases,
                   const std::string& quals);

  const std::string& FragmentId() const { return core_info_.fragment_id; }
  const std::string& Bases() const { return core_info_.bases; }
  const std::string& Quals() const { return core_info_.quals; }

  // Provide access to information fro SAM files.
  int32_t SamChromId() const { return sam_info_.chrom_id; }
  void SetSamChromId(int32_t chrom_id) { sam_info_.chrom_id = chrom_id; }

  int32_t SamPos() const { return sam_info_.pos; }
  void SetSamPos(int32_t pos) { sam_info_.pos = pos; }

  int32_t SamMapq() const { return sam_info_.mapq; }
  void SetSamMapq(int32_t mapq) { sam_info_.mapq = mapq; }

  int32_t SamMateChromId() const { return sam_info_.mate_chrom_id; }
  void SetSamMateChromId(int32_t mate_chrom_id) {
    sam_info_.mate_chrom_id = mate_chrom_id;
  }

  int32_t SamMatePos() const { return sam_info_.mate_pos; }
  void SetSamMatePos(int32_t mate_pos) { sam_info_.mate_pos = mate_pos; }

  bool IsSamMapped() const { return sam_info_.is_mapped; }
  void SetIsSamMapped(bool is_mapped) { sam_info_.is_mapped = is_mapped; }

  bool IsFirstMate() const { return sam_info_.is_first_mate; }
  void SetIsFirstMate(bool is_first_mate) {
    sam_info_.is_first_mate = is_first_mate;
  }

  bool IsMateSamMapped() const { return sam_info_.is_mate_mapped; }
  void SetIsMateSamMapped(bool is_mate_mapped) {
    sam_info_.is_mate_mapped = is_mate_mapped;
  }

  // Provide access to graph-specific information.
  const GraphMapping& CanonicalMapping() const;

  bool HasCanonicalMapping() const {
    return ((bool)graph_info_.canonical_mapping_ptr);
  }

  void SetCanonicalMapping(const GraphMapping& graph_mapping) {
    graph_info_.canonical_mapping_ptr.reset(new GraphMapping(graph_mapping));
  }

  MappingType CanonicalMappingType() const {
    return graph_info_.canonical_mapping_type;
  }

  void SetCanonicalMappingType(MappingType mapping_type) {
    graph_info_.canonical_mapping_type = mapping_type;
  }

  int32_t NumStrUnitsSpanned() const {
    return graph_info_.num_str_units_spanned;
  }

  void SetNumStrUnitsSpanned(int32_t num_str_units_spanned) {
    graph_info_.num_str_units_spanned = num_str_units_spanned;
  }

  bool operator==(const Read& other) const {
    return (core_info_ == other.core_info_ && sam_info_ == other.sam_info_ &&
            graph_info_ == other.graph_info_);
  }

 private:
  CoreInfo core_info_;
  SamInfo sam_info_;
  GraphInfo graph_info_;
};

typedef std::shared_ptr<Read> ReadPtr;

std::ostream& operator<<(std::ostream& os, const Read& read);

}  // namespace reads