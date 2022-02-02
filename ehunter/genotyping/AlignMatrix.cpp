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

#include "genotyping/AlignMatrix.hh"

#include <algorithm>
#include <cassert>

#include <numeric>

using graphtools::GraphAlignment;
using std::vector;

namespace ehunter
{
namespace strgt
{

AlignMatrix::AlignMatrix(int strNode)
    : strNode_(strNode)
    , alignmentCalculator_(strNode)
{
}

StrAlign AlignMatrix::getAlign(int readIndex, int alleleSize) const
{
    if (readIndex >= static_cast<int>(alignScoreMatrix_.size()))
    {
        throw std::runtime_error("Encountered invalid alignment matrix index " + std::to_string(readIndex));
    }

    if (alleleSize < static_cast<int>(alignScoreMatrix_[readIndex].size()))
    {
        return alignScoreMatrix_[readIndex][alleleSize];
    }
    else
    {
        assert(!alignScoreMatrix_[readIndex].empty());
        return alignScoreMatrix_[readIndex].back();
    }
}

void AlignMatrix::add(const GraphAlignment& read, const GraphAlignment& mate)
{
    const int numMotifsInRead = std::count(read.path().begin(), read.path().end(), strNode_);
    const int numMorifsInMate = std::count(mate.path().begin(), mate.path().end(), strNode_);

    if (numMotifsInRead != 0 || numMorifsInMate != 0)
    {
        add(read);
        add(mate);
    }
}

void AlignMatrix::add(const graphtools::GraphAlignment& graphAlign)
{
    vector<StrAlign> strAligns;
    const int numMotifsInAlign = std::count(graphAlign.path().begin(), graphAlign.path().end(), strNode_);
    StrAlign alignToMostConsistentAllele = alignmentCalculator_.findConsistentAlignment(numMotifsInAlign, graphAlign);

    bestAlignsByRead_.push_back(alignToMostConsistentAllele);

    for (int numMotifs = numMotifsInAlign - 1; numMotifs != -1; --numMotifs)
    {
        StrAlign align = alignmentCalculator_.findConsistentAlignment(numMotifs, graphAlign);
        strAligns.emplace_back(align);
    }
    std::reverse(strAligns.begin(), strAligns.end());
    strAligns.emplace_back(alignToMostConsistentAllele);
    StrAlign previousAlign = strAligns.back();
    for (int numMotifs = numMotifsInAlign + 1;; ++numMotifs)
    {
        StrAlign align = alignmentCalculator_.findConsistentAlignment(numMotifs, graphAlign);
        if (align.type() == previousAlign.type() && align.score() == previousAlign.score())
        {
            break;
        }
        strAligns.emplace_back(align);
        previousAlign = align;
    }

    alignScoreMatrix_.emplace_back(strAligns);
}

int AlignMatrix::getMaxMotifCount() const
{
    int maxMotifCount = 0;
    for (const auto& readAligns : alignScoreMatrix_)
    {
        const int maxMotifsInAlign = static_cast<int>(readAligns.size()) - 1;
        if (maxMotifsInAlign > maxMotifCount)
        {
            maxMotifCount = maxMotifsInAlign;
        }
    }

    return maxMotifCount;
}

StrAlign AlignMatrix::getBestAlign(int readIndex) const { return bestAlignsByRead_[readIndex]; }

std::ostream& operator<<(std::ostream& out, const AlignMatrix& matrix)
{
    const auto& alignsByRead = matrix.matrix();

    const unsigned alignCount(alignsByRead.size());
    std::vector<unsigned> alignIndexes(alignCount);
    std::iota(alignIndexes.begin(), alignIndexes.end(), 0);

    // Sort alignIndexes to canonicalize the alignmatrix output
    std::sort(alignIndexes.begin(), alignIndexes.end(), [&alignsByRead](const unsigned ai, const unsigned bi) -> bool {
        const auto& a(alignsByRead[ai]);
        const auto& b(alignsByRead[bi]);
        if (a.size() < b.size())
        {
            return true;
        }
        if (b.size() < a.size())
        {
            return false;
        }
        for (unsigned i(0); i < a.size(); ++i)
        {
            if (a[i] < b[i])
            {
                return true;
            }
            if (b[i] < a[i])
            {
                return false;
            }
        }
        return false;
    });

    auto dumpStrAlign = [&out](const StrAlign& align) {
        out << "(";
        switch (align.type())
        {
        case StrAlign::Type::kOutside:
            out << "O";
            break;
        case StrAlign::Type::kInRepeat:
            out << "I";
            break;
        case StrAlign::Type::kSpanning:
            out << "S";
            break;
        case StrAlign::Type::kFlanking:
            out << "F";
            break;
        }
        out << "," << align.numMotifs() << "," << align.score() << "), ";
    };

    for (unsigned origAlignIndex(0); origAlignIndex < alignCount; ++origAlignIndex)
    {
        const unsigned alignIndex(alignIndexes[origAlignIndex]);
        const auto& aligns(alignsByRead[alignIndex]);
        for (const StrAlign& align : aligns)
        {
            dumpStrAlign(align);
        }
        out << "\n";
    }

    return out;
}

void addIrrPairsIfPossibleExpansion(int maxMotifsInRead, AlignMatrix& alignMatrix, int numIrrPairs)
{
    // Find highest-scoring in-repeat read
    assert(alignMatrix.bestAlignsByRead_.size() == alignMatrix.alignScoreMatrix_.size());

    const int longIrrLowerBound = static_cast<int>(0.90 * maxMotifsInRead);

    int topIrrIndex = -1;
    int topIrrScore = -1;
    for (int readIndex = 0; readIndex != alignMatrix.numReads(); ++readIndex)
    {
        const StrAlign& align = alignMatrix.bestAlignsByRead_[readIndex];
        const bool isLongIrr = align.type() == StrAlign::Type::kInRepeat && align.numMotifs() >= longIrrLowerBound;

        if (isLongIrr && align.score() > topIrrScore)
        {
            topIrrScore = align.score();
            topIrrIndex = readIndex;
        }
    }

    if (topIrrIndex == -1)
    {
        return;
    }

    StrAlign irrTopAlign = alignMatrix.bestAlignsByRead_[topIrrIndex];
    vector<StrAlign> irrAligns = alignMatrix.alignScoreMatrix_[topIrrIndex];

    for (int newIrrIndex = 0; newIrrIndex != 2 * numIrrPairs; ++newIrrIndex)
    {
        alignMatrix.bestAlignsByRead_.push_back(irrTopAlign);
        alignMatrix.alignScoreMatrix_.push_back(irrAligns);
    }
}

void AlignMatrix::remove(int readIndex)
{
    assert(readIndex < numReads());
    bestAlignsByRead_.erase(bestAlignsByRead_.begin() + readIndex);
    alignScoreMatrix_.erase(alignScoreMatrix_.begin() + readIndex);
}

}
}
