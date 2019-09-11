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

#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "thirdparty/intervaltree/IntervalTree.h"

#include "reads/Read.hh"
#include "region/Region.hh"
#include "sample_analysis/GenomeMask.hh"

namespace ehunter
{

// Specifies which mates should be processed with a given locus analyzer
enum class AnalyzerInputType
{
    kReadOnly,
    kMateOnly,
    kBothReads
};

// Stores information needed to properly pass reads to the analyzer
struct AnalyzerBundle
{
    explicit AnalyzerBundle(Region* regionPtr)
        : inputType(AnalyzerInputType::kBothReads)
        , regionPtr(regionPtr)
    {
    }

    AnalyzerInputType inputType;
    Region* regionPtr;
};

// Enables retrieval of appropriate locus analyzers by genomic coordinates of read alignments
class AnalyzerFinder
{
public:
    explicit AnalyzerFinder(std::vector<Region::SPtr>& locusAnalyzerPtrs);

    // Retrieves analyzers appropriate for the given read pair
    std::vector<AnalyzerBundle> query(
        int32_t readContigId, int64_t readStart, int64_t readEnd, int32_t mateContigId, int64_t mateStart,
        int64_t mateEnd) const;

    // Retrieves analyzers appropriate for the given read
    std::vector<AnalyzerBundle> query(int32_t readContigId, int64_t readStart, int64_t readEnd) const;

private:
    using AnalyzerIntervalTree = ehunter::IntervalTree<std::size_t, AnalyzerBundle>;
    using AnalyzerIntervalTrees = std::unordered_map<int32_t, AnalyzerIntervalTree>;

    AnalyzerIntervalTrees intervalTrees_;
};

}
