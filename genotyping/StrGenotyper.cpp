//
// ExpansionHunter
// Copyright 2016-2020 Illumina, Inc.
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

#include "genotyping/StrGenotyper.hh"

#include <vector>

#include "genotyping/AlignMatrixFiltering.hh"
#include "genotyping/OneAlleleStrGenotyper.hh"
#include "genotyping/TwoAlleleStrGenotyper.hh"

using std::unordered_set;
using std::vector;

namespace ehunter
{
namespace strgt
{

unordered_set<int> getAlleleCandidates(int readLen, int motifLen, const AlignMatrix& alignMatrix)
{
    unordered_set<int> candidateSizes;

    int numInRepeatReads = 0;
    int numFlankingReads = 0;
    int longestFlankingSize = 0;

    for (int readIndex = 0; readIndex != static_cast<int>(alignMatrix.numReads()); ++readIndex)
    {
        auto topAlign = alignMatrix.getBestAlign(readIndex);
        if (topAlign.type() == StrAlign::Type::kSpanning)
        {
            candidateSizes.insert(topAlign.numMotifs());
            numFlankingReads += 2;
        }
        else if (topAlign.type() == StrAlign::Type::kFlanking)
        {
            longestFlankingSize = std::max(longestFlankingSize, topAlign.numMotifs());
            ++numFlankingReads;
        }
        else if (topAlign.type() == StrAlign::Type::kInRepeat)
        {
            ++numInRepeatReads;
        }
    }

    if (candidateSizes.empty() || *std::max_element(candidateSizes.begin(), candidateSizes.end()) < longestFlankingSize)
    {
        candidateSizes.insert(longestFlankingSize);
    }

    if (numFlankingReads > 0 && numInRepeatReads > 0)
    {
        candidateSizes.insert(static_cast<int>(static_cast<double>(readLen) / motifLen));
        double depth = static_cast<double>(numFlankingReads) / 2;
        double mediumExpansion = readLen + static_cast<double>(numInRepeatReads * readLen) / depth;
        candidateSizes.insert(static_cast<int>(mediumExpansion / motifLen));
        double longExpansion = readLen + static_cast<double>(2 * numInRepeatReads * readLen) / depth;
        candidateSizes.insert(static_cast<int>(longExpansion / motifLen));
    }

    return candidateSizes;
}

vector<double> getTopFragLogliks(FragLogliks& loglikCalc, const unordered_set<int>& alleleCandidates)
{
    const double negInf = std::numeric_limits<double>::lowest();
    vector<double> topFragLogliks(loglikCalc.numFrags(), negInf);

    for (int fragIndex = 0; fragIndex != loglikCalc.numFrags(); ++fragIndex)
    {
        for (const int alleleSize : alleleCandidates)
        {
            const double fragLoglik = loglikCalc.getLoglik(fragIndex, alleleSize);
            if (fragLoglik > topFragLogliks[fragIndex])
            {
                topFragLogliks[fragIndex] = fragLoglik;
            }
        }
    }

    return topFragLogliks;
}

RepeatGenotype genotype(AlleleCount alleleCount, int motifLen, int readLen, int fragLen, AlignMatrix& alignMatrix)
{
    filter(alignMatrix);
    FragLogliks fragLoglikCalc(motifLen, readLen, fragLen, &alignMatrix);
    unordered_set<int> candidateAlleleSizes = getAlleleCandidates(readLen, motifLen, alignMatrix);
    vector<double> topFragLogliks = getTopFragLogliks(fragLoglikCalc, candidateAlleleSizes);

    if (alleleCount == AlleleCount::kTwo)
    {
        TwoAlleleGenotyper genotyper(motifLen, fragLen, std::move(topFragLogliks), &fragLoglikCalc);
        return genotyper.genotype(candidateAlleleSizes);
    }
    else
    {
        assert(alleleCount == AlleleCount::kOne);
        OneAlleleGenotyper genotyper(motifLen, std::move(topFragLogliks), &fragLoglikCalc);
        return genotyper.genotype(candidateAlleleSizes);
    }
}

}
}