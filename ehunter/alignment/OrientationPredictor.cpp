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

#include "alignment/OrientationPredictor.hh"

#include <algorithm>
#include <iostream>
#include <list>
#include <vector>

#include "graphcore/Path.hh"
#include "graphcore/PathOperations.hh"
#include "graphutils/SequenceOperations.hh"

using graphtools::KmerIndex;
using graphtools::Path;
using std::list;
using std::string;
using std::vector;

namespace ehunter
{

static int countNonoverlappingKmerMatches(const string& query, const KmerIndex& kmerIndex)
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
        = countNonoverlappingKmerMatches(graphtools::reverseComplement(query), kmerIndex_);

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
        out << "kAlignsInOriginalOrientation";
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
