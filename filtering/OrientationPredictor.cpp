//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
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

#include "filtering/OrientationPredictor.hh"

#include <algorithm>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "graphalign/GaplessAligner.hh"
#include "graphcore/Path.hh"
#include "graphcore/PathOperations.hh"

using graphtools::alignWithoutGaps;
using graphtools::GraphAlignment;
using graphtools::KmerIndex;
using graphtools::Path;
using std::list;
using std::string;
using std::vector;

namespace ehunter
{

static int countNonoverlappingKmerMatches(const string query, const KmerIndex& kmerIndex)
{
    const std::size_t kmerLength = kmerIndex.kmerLength();
    int matchCount = 0;
    size_t position = 0;
    while (position + kmerLength <= query.length())
    {
        string kmer = query.substr(position, kmerLength);
        std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);

        if (kmerIndex.contains(kmer))
        {
            ++matchCount;
            position += kmerLength;
        }
        else
        {
            ++position;
        }
    }
    return matchCount;
}

OrientationPrediction OrientationPredictor::predict(const std::string& query) const
{
    const int numForwardMatches = countNonoverlappingKmerMatches(query, kmerIndex_);
    const int numReverseComplementMatches
        = countNonoverlappingKmerMatches(query, kmerIndexForReverseComplementedGraph_);

    const int maxMatches = std::max(numForwardMatches, numReverseComplementMatches);

    if (maxMatches < minKmerMatchesToPass_)
    {
        return OrientationPrediction::kDoesNotAlign;
    }

    if (numForwardMatches >= numReverseComplementMatches)
    {
        return OrientationPrediction::kAlignsInOriginalOrientation;
    }
    else
    {
        return OrientationPrediction::kAlignsInReverseComplementOrientation;
    }
}

std::ostream& operator<<(std::ostream& out, OrientationPrediction orientationPrediction)
{
    switch (orientationPrediction)
    {
    case OrientationPrediction::kAlignsInOriginalOrientation:
        out << "kAlignsInReverseComplementOrientation";
        break;
    case OrientationPrediction::kAlignsInReverseComplementOrientation:
        out << "kAlignsInReverseComplementOrientation";
        break;
    case OrientationPrediction::kDoesNotAlign:
        out << "kDoesNotAlign";
    }

    return out;
}

}
