//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "reads/read.h"

#include <stdexcept>

using std::string;

namespace reads {
void Read::SetCoreInfo(const string& fragment_id, const string& bases,
                       const string& quals) {
  if (bases.length() != quals.length()) {
    throw std::runtime_error("Bases " + bases + " and quals " + quals +
                             " must have the same length");
  }

  core_info_.fragment_id = fragment_id;
  core_info_.bases = bases;
  core_info_.quals = quals;
}

const string& Read::FragmentId() const {
  if (core_info_.fragment_id.empty()) {
    throw std::logic_error("Cannot access fragment id that was not set");
  }
  return core_info_.fragment_id;
}
const string& Read::Bases() const {
  if (core_info_.fragment_id.empty()) {
    throw std::logic_error("Cannot access read sequence that was not set");
  }
  return core_info_.bases;
}
const string& Read::Quals() const {
  if (core_info_.fragment_id.empty()) {
    throw std::logic_error("Cannot access read quals that was not set");
  }
  return core_info_.quals;
}

const GraphMapping& Read::CanonicalMapping() const {
  if (!graph_info_.canonical_mapping_ptr) {
    throw std::logic_error(
        "Attempted access to uninitialized canonical mapping");
  }
  return *graph_info_.canonical_mapping_ptr;
}

}  // namespace reads