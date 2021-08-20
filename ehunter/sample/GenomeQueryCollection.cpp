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

#include "sample/GenomeQueryCollection.hh"

using ehunter::locus::LocusAnalyzer;
using std::unique_ptr;
using std::vector;

namespace ehunter
{

namespace
{
void initializeGenomeMask(GenomeMask& genomeMask, vector<unique_ptr<locus::LocusAnalyzer>>& locusAnalyzers)
{
    for (auto& locusAnalyzer : locusAnalyzers)
    {
        const LocusSpecification& locusSpec = locusAnalyzer->locusSpec();

        for (const auto& region : locusSpec.targetReadExtractionRegions())
        {
            genomeMask.addRegion(region.contigIndex(), region.start(), region.end());
        }

        for (const auto& region : locusSpec.offtargetReadExtractionRegions())
        {
            genomeMask.addRegion(region.contigIndex(), region.start(), region.end());
        }
    }
}
}

GenomeQueryCollection::GenomeQueryCollection(vector<unique_ptr<LocusAnalyzer>>& locusAnalyzers)
    : analyzerFinder(locusAnalyzers)
{
    initializeGenomeMask(targetRegionMask, locusAnalyzers);
}

}
