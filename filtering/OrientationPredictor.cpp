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

#include "filtering/OrientationPredictor.hh"

#include <algorithm>
#include <iostream>
#include <list>
#include <stdexcept>
#include <string>
#include <vector>
#include <stdexcept>

#include "graphcore/Path.hh"

using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::Path;
using graphtools::reverseGraph;
using std::list;
using std::string;
using std::vector;

namespace ehunter
{

static int countNonoverlappingKmerMatches(const string& query, int kmerLength, const BloomFilter& bloomFilter)
{
    int matchCount = 0;
    size_t position = 0;
    while (position + kmerLength <= query.length())
    {
        string kmer = query.substr(position, kmerLength);
        std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);

        if (bloomFilter.maybeContains(kmer))
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

OrientationPredictor::OrientationPredictor(const Graph* graph)
    : kmerLength_(10)
    , minKmerMatchesToPass_(3)
    , bloomFilter_(build(*graph, kmerLength_))
    , oppositeBloomFilter_(build(reverseGraph(*graph, true), kmerLength_))

{
}

OrientationPrediction OrientationPredictor::predict(const string& query) const
{
    const int numMatches = countNonoverlappingKmerMatches(query, kmerLength_, bloomFilter_);
    const int numOppositeMatches = countNonoverlappingKmerMatches(query, kmerLength_, oppositeBloomFilter_);

    const int maxMatches = std::max(numMatches, numOppositeMatches);

    if (maxMatches < minKmerMatchesToPass_)
    {
        return OrientationPrediction::kDoesNotAlign;
    }

    if (numMatches >= numOppositeMatches)
    {
        return OrientationPrediction::kAlignsInOriginalOrientation;
    }
    else
    {
        return OrientationPrediction::kAlignsInOppositeOrientation;
    }
}

std::ostream& operator<<(std::ostream& out, OrientationPrediction orientationPrediction)
{
    switch (orientationPrediction)
    {
    case OrientationPrediction::kAlignsInOriginalOrientation:
        out << "kAlignsInOriginalOrientation";
        break;
    case OrientationPrediction::kAlignsInOppositeOrientation:
        out << "kAlignsInOppositeOrientation";
        break;
    case OrientationPrediction::kDoesNotAlign:
        out << "kDoesNotAlign";
    }

    return out;
}

}
