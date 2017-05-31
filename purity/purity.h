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

#include <vector>
#include <string>

std::vector<std::vector<std::string>> shift_units(
    const std::vector<std::string>& units);

double MatchRepeatRc(const std::vector<std::vector<std::string>>& units_shifts,
                     const std::string& bases, const std::string& quals,
                     size_t min_baseq = 20);

double MatchRepeat(const std::vector<std::vector<std::string>>& units_shifts,
                   const std::string& bases, const std::string& quals,
                   size_t& match_offset, size_t min_baseq = 20);

double MatchRepeat(const std::vector<std::string>& units,
                   const std::string& bases, const std::string& quals,
                   size_t min_baseq = 20);

double MatchUnits(const std::vector<std::string>& units,
                  std::string::const_iterator bases_start,
                  std::string::const_iterator bases_end,
                  std::string::const_iterator quals_start,
                  std::string::const_iterator quals_end, size_t min_baseq = 20);
