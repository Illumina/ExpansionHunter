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
#include <string>

#include "graphs/graph.h"
#include "graphs/graph_mapping.h"

void SplitNodeCigar(const std::string& node_cigar, std::string& cigar,
                    int32_t& node_id);

GraphMapping DecodeFromString(int32_t first_node_start,
                              const std::string& graph_cigar,
                              const std::string& query,
                              GraphSharedPtr graph_ptr);

std::string EncodeGraphMapping(const GraphMapping& graph_mapping,
                               int32_t padding = 0);