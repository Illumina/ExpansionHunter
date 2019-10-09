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
#include "sample_analysis/GenomeMask.hh"
#include "workflow/RegionModel.hh"

namespace ehunter
{

// Enables retrieval of appropriate locus analyzers by genomic coordinates of read alignments
class ModelFinder
{
public:
    explicit ModelFinder(const std::vector<std::shared_ptr<RegionModel>>& models);

    // Retrieves models whose regions overlap the given interval
    std::unordered_set<RegionModel*> query(int contigId, int64_t start, int64_t end) const;

private:
    using ModelTree = ehunter::IntervalTree<std::size_t, RegionModel*>;
    using ContigToModelTree = std::unordered_map<int, ModelTree>;

    ContigToModelTree contigToModelTrees_;
};

}
