//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Felix Schlesinger <fschlesinger@illumina.com>
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

#include <unordered_set>

#include "reads/Read.hh"
#include "region/RegionModel.hh"

namespace ehunter
{

void dispatch(
    int readContig, int64_t readStart, int64_t readEnd, const Read& read, int mateContig, int64_t mateStart,
    int64_t mateEnd, const Read& mate, const std::unordered_set<RegionModel*>& models);

}