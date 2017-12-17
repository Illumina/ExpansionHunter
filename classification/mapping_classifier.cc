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

#include "classification/mapping_classifier.h"

using std::map;
using std::ostream;
using std::string;

ostream& operator<<(ostream& os, const ReadClass& read_class) {
  static const std::map<ReadClass, string> class_to_string = {
      {ReadClass::kSpansRepeat, "kSpansRepeat"},
      {ReadClass::kFlanksRepeat, "kFlanksRepeat"},
      {ReadClass::kInsideRepeat, "kInsideRepeat"},
      {ReadClass::kOutsideRepeat, "kOutsideRepeat"},
      {ReadClass::kUnmapped, "kUnmapped"},
      {ReadClass::kUnknown, "kUnknown"}};
  os << class_to_string.at(read_class);
  return os;
}

ReadClass StrMappingClassifier::Classify(const GraphMapping& mapping) {
  const bool overlaps_left_flank = mapping.OverlapsNode(left_flank_id_);
  const bool overlaps_repeat_unit = mapping.OverlapsNode(repeat_unit_id_);
  const bool overlaps_right_flank = mapping.OverlapsNode(right_flank_id_);
  const bool overlaps_both_flanks = overlaps_left_flank && overlaps_right_flank;
  const bool overlaps_either_flank =
      overlaps_left_flank || overlaps_right_flank;

  if (overlaps_both_flanks) {
    return ReadClass::kSpansRepeat;
  }

  if (overlaps_either_flank && overlaps_repeat_unit) {
    return ReadClass::kFlanksRepeat;
  }

  if (overlaps_repeat_unit) {
    return ReadClass::kInsideRepeat;
  }

  if (overlaps_either_flank) {
    return ReadClass::kOutsideRepeat;
  }

  return ReadClass::kUnmapped;
}
