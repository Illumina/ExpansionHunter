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

#pragma once

#include <map>
#include <string>
#include <vector>

#include "common/common.h"
#include "common/parameters.h"
#include "common/repeat_spec.h"
#include "genotyping/short_repeat_genotyper.h"

void GenotypeRepeat(const Parameters &parameters, const RepeatSpec &repeat_spec,
                    int max_num_units_in_read, double prop_correct_molecules,
                    double hap_depth, int read_len,
                    const std::vector<RepeatAllele> &haplotype_candidates,
                    const std::map<int, int> &flanking_size_count,
                    const std::map<int, int> &spanning_size_count,
                    GenotypeType genotype_type, RepeatGenotype &genotype);