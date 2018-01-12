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

#include "classification/overlap_quantification.h"

#include <list>

using std::list;

int32_t StrOverlapQuantifier::NumUnitsOverlapped(
    const GraphMapping& mapping) const {
  const list<int32_t> repeat_unit_indexes =
      mapping.GetIndexesOfNode(repeat_unit_id_);
  int32_t num_units_overlapped =
      static_cast<int32_t>(repeat_unit_indexes.size());
  if (!repeat_unit_indexes.empty()) {
    const int32_t last_index = repeat_unit_indexes.back();
    const Mapping& last_mapping = mapping[last_index];
    const int32_t alignment_span = last_mapping.ReferenceSpan();

    if (str_unit_len_ != alignment_span) {
      --num_units_overlapped;
    }
  }

  return num_units_overlapped;
}
