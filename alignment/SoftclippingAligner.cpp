//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "alignment/SoftclippingAligner.hh"

#include "alignment/GraphAlignmentOperations.hh"
#include "alignment/HighQualityBaseRunFinder.hh"

using graphtools::Graph;
using graphtools::GraphAlignment;
using std::list;
using std::string;

namespace ehunter
{

SoftclippingAligner::SoftclippingAligner(
    const Graph* graphPtr, const std::string& alignerName, int kmerLenForAlignment, int paddingLength,
    int seedAffixTrimLength)
    : aligner_(graphPtr, kmerLenForAlignment, paddingLength, seedAffixTrimLength, alignerName)
{
}

list<GraphAlignment> SoftclippingAligner::align(const string& query) const
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

    return aligner_.align(query);
}

}
