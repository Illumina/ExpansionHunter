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

#include "alignment/SoftclippingAligner.hh"

#include "alignment/HighQualityBaseRunFinder.hh"
#include "alignment/OperationsOnAlignments.hh"

using graphtools::Graph;
using graphtools::GraphAlignment;
using std::list;
using std::string;

namespace ehunter
{

SoftclippingAligner::SoftclippingAligner(
    const Graph* graphPtr, int kmerLenForAlignment, int paddingLength, int seedAffixTrimLength)
    : aligner_(graphPtr, kmerLenForAlignment, paddingLength, seedAffixTrimLength)
{
}

list<GraphAlignment> SoftclippingAligner::align(const string& query, graphtools::AlignerSelector& alignerSelector) const
{
    /*
    HighQualityBaseRunFinder goodBaseFinder(6, 2, query.size() / 2);
    const auto goodBasesRange = goodBaseFinder.find(query);
    const auto goodBases = string(goodBasesRange.first, goodBasesRange.second);

    int numBasesTrimmedFromLeft = goodBasesRange.first - query.begin();
    int numBasesTrimmedFromRight = query.end() - goodBasesRange.second;

    list<GraphAlignment> extendedAlignments;
    for (const auto& alignment : aligner_.align(goodBases))
    {
        const auto extendedAlignment = extendWithSoftclip(alignment, numBasesTrimmedFromLeft, numBasesTrimmedFromRight);
        extendedAlignments.push_back(extendedAlignment);
    }
    */

    return aligner_.align(query, alignerSelector);
}

}
