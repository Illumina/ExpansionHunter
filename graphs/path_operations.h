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

#include <list>
#include <string>
#include <vector>

#include "graphs/path.h"

/**
 * @brief Splits sequence into segments corresponding to the path
 *
 * @param path Any valid path
 * @param sequence A string having the same length as the path
 * @return Segments of the sequence corresponding to nodes spanned by the path
 */
std::vector<std::string> SplitByPath(const GraphPath& path,
                                     const std::string& sequence);

std::list<GraphPath> ComputeRightEndings(const GraphPath& path,
                                         int32_t dist_from_right_end);

std::list<GraphPath> ComputeLeftEndings(const GraphPath& path,
                                        int32_t dist_from_left_end);
