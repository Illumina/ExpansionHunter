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

#include <string>

#include "graphs/graph.h"

Graph makeDeletionGraph(const std::string& left_flank,
                        const std::string& deletion,
                        const std::string& right_flank);

Graph makeSwapGraph(const std::string& left_flank, const std::string& deletion,
                    const std::string& insertion,
                    const std::string& right_flank);

Graph makeDoubleSwapGraph(const std::string& left_flank,
                          const std::string& deletion1,
                          const std::string& insertion1,
                          const std::string& middle,
                          const std::string& deletion2,
                          const std::string& insertion2,
                          const std::string& right_flank);

Graph makeLooplessStrGraph(int32_t read_len, const std::string& left_flank,
                           const std::string& repeat_unit,
                           const std::string& right_flank);
