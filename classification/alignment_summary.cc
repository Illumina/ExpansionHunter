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

#include "classification/alignment_summary.h"

#include <stdexcept>

using std::map;
using std::vector;

void SummarizeAlignments(const vector<reads::ReadPtr> &read_ptrs,
                         map<int32_t, int32_t> &flanking_size_counts,
                         map<int32_t, int32_t> &spanning_size_counts) {
  for (const reads::ReadPtr &read_ptr : read_ptrs) {
    const MappingType read_type = read_ptr->CanonicalMappingType();
    const int32_t num_str_units_spanned = read_ptr->NumStrUnitsSpanned();
    switch (read_type) {
      case MappingType::kFlanksRepeat:
      case MappingType::kInsideRepeat:
        flanking_size_counts[num_str_units_spanned] += 1;
        break;
      case MappingType::kSpansRepeat:
        spanning_size_counts[num_str_units_spanned] += 1;
      default:
        break;
    }
  }
}