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

#include <memory>
#include <string>
#include <unordered_map>

#include "core/LocusStats.hh"
#include "locus/VariantFindings.hh"

namespace ehunter
{

// Container with per-locus analysis results
struct LocusFindings
{
    explicit LocusFindings(LocusStats stats = {})
        : stats(stats)
    {
    }
    LocusStats stats;
    // VariantFindings is an abstract class from which findings for all variant types are derived
    std::unordered_map<std::string, std::unique_ptr<VariantFindings>> findingsForEachVariant;
};

using SampleFindings = std::vector<LocusFindings>;

}
