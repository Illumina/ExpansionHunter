//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
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
