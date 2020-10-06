//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>
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
#include <vector>

namespace ehunter
{

/// Binned Bit-mask covering the genome
/// Used to mark target regions for fast read screening
class GenomeMask
{
public:
    void addRegion(int32_t contigId, int64_t start, int64_t stop);
    bool query(int32_t contigId, int64_t pos) const;

private:
    using contigMask = std::vector<bool>;
    std::vector<contigMask> mask_;
};

}
