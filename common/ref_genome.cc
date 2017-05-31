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

#include "common/ref_genome.h"

#include <algorithm>
#include <string>
#include <stdexcept>

using std::string;

RefGenome::RefGenome(const string& genome_path) : genome_path_(genome_path) {
  fai_ptr_ = fai_load(genome_path_.c_str());
}

RefGenome::~RefGenome() { fai_destroy(fai_ptr_); }

// Load reference sequence specified by region.
void RefGenome::ExtractSeq(const string& region, string* sequence) const {
  int len;  // throwaway...

  char* ref_tmp = fai_fetch(fai_ptr_, region.c_str(), &len);

  if (!ref_tmp || len == -1 || len == -2) {
    throw std::runtime_error("ERROR: can't extract " + region + " from "
        + genome_path_ + "; in particular, chromosome names must match "
        "exactly (e.g. \"chr1\" and \"1\" are distinct names)");
  }

  sequence->assign(ref_tmp);
  free(ref_tmp);

  std::transform(sequence->begin(), sequence->end(), sequence->begin(),
                 ::toupper);
}