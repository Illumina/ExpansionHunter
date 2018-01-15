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

#include <map>
#include <string>
#include <unordered_set>

#include "classification/mapping_classifier.h"

using std::list;
using std::ostream;
using std::string;
using std::unordered_set;

ostream& operator<<(ostream& os, const MappingType& read_class) {
  static const std::map<MappingType, string> class_to_string = {
      {MappingType::kSpansRepeat, "kSpansRepeat"},
      {MappingType::kFlanksRepeat, "kFlanksRepeat"},
      {MappingType::kInsideRepeat, "kInsideRepeat"},
      {MappingType::kOutsideRepeat, "kOutsideRepeat"},
      {MappingType::kUnmapped, "kUnmapped"}};
  os << class_to_string.at(read_class);
  return os;
}

GraphMapping StrMappingClassifier::GetCanonicalMapping(
    const list<GraphMapping>& mappings) const {
  const GraphMapping* canonical_mapping_ptr = nullptr;
  for (const GraphMapping& mapping : mappings) {
    MappingType mapping_type = Classify(mapping);
    if (!canonical_mapping_ptr) {
      canonical_mapping_ptr = &mapping;
    } else if (mapping_type == MappingType::kInsideRepeat) {
      return mapping;
    } else if (mapping_type == MappingType::kFlanksRepeat) {
      canonical_mapping_ptr = &mapping;
    }
  }
  return *canonical_mapping_ptr;
}

MappingType StrMappingClassifier::Classify(const GraphMapping& mapping) const {
  const bool overlaps_left_flank = mapping.OverlapsNode(left_flank_id_);
  const bool overlaps_repeat_unit = mapping.OverlapsNode(repeat_unit_id_);
  const bool overlaps_right_flank = mapping.OverlapsNode(right_flank_id_);
  const bool overlaps_both_flanks = overlaps_left_flank && overlaps_right_flank;
  const bool overlaps_either_flank =
      overlaps_left_flank || overlaps_right_flank;

  if (overlaps_both_flanks) {
    return MappingType::kSpansRepeat;
  }

  if (overlaps_either_flank && overlaps_repeat_unit) {
    return MappingType::kFlanksRepeat;
  }

  if (overlaps_repeat_unit) {
    return MappingType::kInsideRepeat;
  }

  if (overlaps_either_flank) {
    return MappingType::kOutsideRepeat;
  }

  return MappingType::kUnmapped;
}
