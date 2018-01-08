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
#include <vector>

#include "common/genomic_region.h"
#include "common/repeat_spec.h"
#include "graphs/graph.h"
#include "reads/read_pairs.h"
#include "reads/read_reader.h"

void ExtractReads(const Region& target_region, reads::ReadReader& read_reader,
                  reads::ReadPairs& read_pairs);

void ExtractReads(const std::vector<Region>& target_regions,
                  reads::ReadReader& read_reader, reads::ReadPairs& read_pairs);

void ExtractReads(const RepeatSpec& repeat_spec, int32_t extension_len,
                  reads::ReadReader& read_reader, reads::ReadPairs& read_pairs);

void AlignReads(const std::shared_ptr<Graph>& graph_ptr, int32_t kmer_len,
                std::vector<reads::ReadPtr>& reads);
