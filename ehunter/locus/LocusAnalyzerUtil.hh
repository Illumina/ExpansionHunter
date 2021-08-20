//
// ExpansionHunter
// Copyright 2016-2021 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//         Chris Saunders <csaunders@illumina.com>
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

#include "locus/LocusAnalyzer.hh"

namespace ehunter
{
namespace locus
{

/// Initialize a LocusAnalyzer for each locus in \p regionCatalog
///
/// \param[in] threadCount Number of threads to distribute initialization over
///
std::vector<std::unique_ptr<LocusAnalyzer>> initializeLocusAnalyzers(
    const RegionCatalog& regionCatalog, const HeuristicParameters& heuristicParams, AlignWriterPtr alignmentWriter,
    int threadCount);

}
}
