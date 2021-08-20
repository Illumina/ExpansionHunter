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

#include "AnalyzerFinder.hh"
#include "GenomeMask.hh"

namespace ehunter
{

// Aggregates various methods for querying genome
struct GenomeQueryCollection
{
    GenomeQueryCollection(std::vector<std::unique_ptr<locus::LocusAnalyzer>>& locusAnalyzers);

    AnalyzerFinder analyzerFinder; // Analyzers searchable by targeted region
    GenomeMask targetRegionMask; // Marks targeted regions to enable fast read screening
};

}
