//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
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
#include <utility>

namespace ehunter
{

using StringIterPair = std::pair<std::string::const_iterator, std::string::const_iterator>;

/**
 * Searches for the first sufficiently-long run of high-quality bases
 *
 * @param query: any query sequence
 * @param probOfGoodBaseInBadRun: probability of observing a high-quality base in a low quality stretch of bases
 * @param probOfGoodBaseInGoodRun: probability of observing a high-quality base in a good quality stretch of bases
 * @return pair of iterators delineating a substring consisting of high quality bases
 */
StringIterPair findHighQualityBaseRun(
    const std::string& query, double probOfGoodBaseInBadRun = 0.1, double probOfGoodBaseInGoodRun = 0.8);

}
