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

#include "genotyping/StrAlign.hh"

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "spdlog/spdlog.h"

#include "graphalign/LinearAlignment.hh"

using graphtools::Alignment;
using graphtools::GraphAlignment;
using graphtools::Operation;
using graphtools::OperationType;

namespace ehunter
{

std::ostream& operator<<(std::ostream& out, StrAlign::Type type)
{
    switch (type)
    {
    case StrAlign::Type::kFlanking:
        out << "StrAlign::Type::kFlanking";
        break;
    case StrAlign::Type::kInRepeat:
        out << "StrAlign::Type::kInRepeat";
        break;
    case StrAlign::Type::kSpanning:
        out << "StrAlign::Type::kSpanning";
        break;
    case StrAlign::Type::kOutside:
        out << "StrAlign::Type::kOutside";
        break;
    default:
        throw std::runtime_error("Encountered unknown StrAlign::Type");
    }

    return out;
}

std::ostream& operator<<(std::ostream& out, const StrAlign& summary)
{
    out << "StrAlign(" << summary.type() << ", " << summary.numMotifs() << ", " << summary.score() << ", "
        << summary.numIndels() << ")";
    return out;
}

void scoreAlignment(
    const Alignment& alignment, int matchScore, int mismatchScore, int gapScore, int& score, int& indelCount)
{
    score = 0;
    indelCount = 0;

    for (const Operation& operation : alignment)
    {
        switch (operation.type())
        {
        case OperationType::kMatch:
            score += matchScore * static_cast<int>(operation.referenceLength());
            break;
        case OperationType::kMismatch:
            score += mismatchScore * static_cast<int>(operation.referenceLength());
            break;
        case OperationType::kInsertionToRef:
            score += gapScore * static_cast<int>(operation.queryLength());
            indelCount += operation.queryLength();
            break;
        case OperationType::kDeletionFromRef:
            score += gapScore * static_cast<int>(operation.referenceLength());
            indelCount += operation.referenceLength();
            break;
        default:
            break;
        }
    }
}

StrAlign ConsistentAlignmentCalculator::clipFromLeft(int numMotifsInAllele, const GraphAlignment& alignment) const
{
    int leftFlankScore = 0, strScore = 0, rightFlankScore = 0;
    int strIndelCount = 0;
    int numMotifsInAlignment = std::count(alignment.path().begin(), alignment.path().end(), strNodeId_);
    int numMotifsLeft = numMotifsInAlignment;
    for (int nodeIndex = 0; nodeIndex != static_cast<int>(alignment.size()); ++nodeIndex)
    {
        int node = alignment.getNodeIdByIndex(nodeIndex);
        const auto& nodeAlign = alignment.alignments()[nodeIndex];
        int nodeScore = 0;
        int nodeIndelCount = 0;
        scoreAlignment(nodeAlign, matchScore_, mismatchScore_, gapOpenScore_, nodeScore, nodeIndelCount);

        if (node < strNodeId_)
        {
            leftFlankScore += nodeScore;
        }
        else if (strNodeId_ < node)
        {
            rightFlankScore += nodeScore;
        }
        else if (numMotifsLeft <= numMotifsInAllele)
        {
            strScore += nodeScore;
            strIndelCount += nodeIndelCount;
            --numMotifsLeft;
        }
        else
        {
            --numMotifsLeft;
        }
    }

    // Zero out negative scores
    leftFlankScore = std::max(leftFlankScore, 0);
    rightFlankScore = std::max(rightFlankScore, 0);

    // Alignment does not overlap the repeat
    if (numMotifsInAlignment == 0 && (leftFlankScore == 0 || rightFlankScore == 0))
    {
        return { StrAlign::Type::kOutside, 0, leftFlankScore + rightFlankScore, 0 };
    }

    const int numCompatibleMotifs = std::min(numMotifsInAlignment, numMotifsInAllele);

    // Original alignment is in-repeat
    if (leftFlankScore == 0 && rightFlankScore == 0)
    {
        return { StrAlign::Type::kInRepeat, numCompatibleMotifs, strScore, strIndelCount };
    }

    // Original alignment is spanning
    if (leftFlankScore > 0 && rightFlankScore > 0)
    {
        if (numMotifsInAlignment == numMotifsInAllele)
        {
            return { StrAlign::Type::kSpanning, numCompatibleMotifs, leftFlankScore + strScore + rightFlankScore,
                     strIndelCount };
        }
        else
        {
            return { StrAlign::Type::kFlanking, numCompatibleMotifs, strScore + rightFlankScore, strIndelCount };
        }
    }

    // Original alignment is right flanking
    if (leftFlankScore == 0 && rightFlankScore > 0)
    {
        return { StrAlign::Type::kFlanking, numCompatibleMotifs, strScore + rightFlankScore, strIndelCount };
    }

    // Original alignment is left flanking
    if (leftFlankScore > 0 && rightFlankScore == 0)
    {
        if (numMotifsInAlignment <= numMotifsInAllele)
        {
            return { StrAlign::Type::kFlanking, numCompatibleMotifs, leftFlankScore + strScore, strIndelCount };
        }
        else
        {
            return { StrAlign::Type::kInRepeat, numCompatibleMotifs, strScore, strIndelCount };
        }
    }

    std::ostringstream out;
    out << "Cannot summarize " << alignment << " clipped from left for STR on node " << strNodeId_;
    spdlog::warn(out.str());

    return { StrAlign::Type::kOutside, 0, 0, 0 };
}

StrAlign
ConsistentAlignmentCalculator::clipFromRight(int numMotifsInAllele, const graphtools::GraphAlignment& alignment) const
{
    int leftFlankScore = 0, strScore = 0, rightFlankScore = 0;
    int strIndelCount = 0;

    int motifIndex = 0;
    for (int nodeIndex = 0; nodeIndex != static_cast<int>(alignment.size()); ++nodeIndex)
    {
        int node = alignment.getNodeIdByIndex(nodeIndex);
        const auto& nodeAlign = alignment.alignments()[nodeIndex];

        int nodeScore = 0;
        int nodeIndelCount = 0;
        scoreAlignment(nodeAlign, matchScore_, mismatchScore_, gapOpenScore_, nodeScore, nodeIndelCount);

        if (node < strNodeId_)
        {
            leftFlankScore += nodeScore;
        }
        else if (strNodeId_ < node)
        {
            rightFlankScore += nodeScore;
        }
        else if (motifIndex + 1 <= numMotifsInAllele)
        {
            assert(node == strNodeId_);
            strScore += nodeScore;
            strIndelCount += nodeIndelCount;
            ++motifIndex;
        }
        else
        {
            assert(node == strNodeId_);
            ++motifIndex;
        }
    }

    // Zero out negative scores
    leftFlankScore = std::max(leftFlankScore, 0);
    rightFlankScore = std::max(rightFlankScore, 0);

    int numMotifsInAlignment = std::count(alignment.path().begin(), alignment.path().end(), strNodeId_);
    // Alignment does not overlap the repeat
    if (numMotifsInAlignment == 0 && (leftFlankScore == 0 || rightFlankScore == 0))
    {
        const int score = leftFlankScore + rightFlankScore;
        return { StrAlign::Type::kOutside, 0, score, 0 };
    }

    const int numCompatibleMotifs = std::min(numMotifsInAlignment, numMotifsInAllele);

    // Original alignment is in-repeat
    if (leftFlankScore == 0 && rightFlankScore == 0)
    {
        return { StrAlign::Type::kInRepeat, numCompatibleMotifs, strScore, strIndelCount };
    }

    // Original alignment is spanning
    if (leftFlankScore > 0 && rightFlankScore > 0)
    {
        if (numMotifsInAlignment == numMotifsInAllele)
        {
            const int score = leftFlankScore + strScore + rightFlankScore;
            return { StrAlign::Type::kSpanning, numCompatibleMotifs, score, strIndelCount };
        }
        else
        {
            const int score = leftFlankScore + strScore;
            return { StrAlign::Type::kFlanking, numCompatibleMotifs, score, strIndelCount };
        }
    }

    // Original alignment is left flanking
    if (leftFlankScore > 0 && rightFlankScore == 0)
    {
        const int score = leftFlankScore + strScore;
        return { StrAlign::Type::kFlanking, numCompatibleMotifs, score, strIndelCount };
    }

    // Original alignment is right flanking
    if (leftFlankScore == 0 && rightFlankScore > 0)
    {
        if (numMotifsInAlignment <= numMotifsInAllele)
        {
            const int score = strScore + rightFlankScore;
            return { StrAlign::Type::kFlanking, numCompatibleMotifs, score, strIndelCount };
        }
        else
        {
            return { StrAlign::Type::kInRepeat, numCompatibleMotifs, strScore, strIndelCount };
        }
    }

    std::ostringstream out;
    out << "Cannot summarize " << alignment << " clipped from right for STR on node " << strNodeId_;
    spdlog::warn(out.str());

    return { StrAlign::Type::kOutside, 0, 0, 0 };
}

StrAlign ConsistentAlignmentCalculator::removeStutter(int numMotifsInAllele, const GraphAlignment& alignment) const
{
    int leftFlankScore = 0, strScore = 0, rightFlankScore = 0;
    int strIndelCount = 0;
    int motifIndex = 0;
    for (int nodeIndex = 0; nodeIndex != static_cast<int>(alignment.size()); ++nodeIndex)
    {
        int node = alignment.getNodeIdByIndex(nodeIndex);
        const auto& nodeAlign = alignment.alignments()[nodeIndex];
        int nodeScore = 0;
        int nodeIndelCount = 0;
        scoreAlignment(nodeAlign, matchScore_, mismatchScore_, gapOpenScore_, nodeScore, nodeIndelCount);

        if (node < strNodeId_)
        {
            leftFlankScore += nodeScore;
        }
        else if (strNodeId_ < node)
        {
            rightFlankScore += nodeScore;
        }
        else if (motifIndex + 1 <= numMotifsInAllele)
        {
            strScore += nodeScore;
            strIndelCount += nodeIndelCount;
            ++motifIndex;
        }
        else
        {
            ++motifIndex;
        }
    }

    // Zero out negative scores
    leftFlankScore = std::max(leftFlankScore, 0);
    rightFlankScore = std::max(rightFlankScore, 0);

    if (leftFlankScore == 0 || rightFlankScore == 0)
    {
        return { StrAlign::Type::kOutside, 0, 0, 0 };
    }

    const int numMotifsInAlignment = std::count(alignment.path().begin(), alignment.path().end(), strNodeId_);
    const int numDiscrepantMotifs = std::abs(numMotifsInAlignment - numMotifsInAllele);
    const int motifLength = alignment.path().graphRawPtr()->nodeSeq(strNodeId_).length();
    const int discrepantLength = motifLength * numDiscrepantMotifs;
    const int gapOpenScore = -24;
    const int gapExtendScore = -12;
    int penaltyScore = numDiscrepantMotifs > 0 ? gapOpenScore + gapExtendScore * (discrepantLength - 1) : 0;
    const int alignmentScore = std::max(leftFlankScore + strScore + penaltyScore + rightFlankScore, 0);

    return { StrAlign::Type::kSpanning, numMotifsInAllele, alignmentScore, strIndelCount };
}

StrAlign
ConsistentAlignmentCalculator::findConsistentAlignment(int numMotifsInAllele, const GraphAlignment& alignment) const
{
    StrAlign stutterFreeAlign = removeStutter(numMotifsInAllele, alignment);
    StrAlign leftClipAlign = clipFromLeft(numMotifsInAllele, alignment);
    StrAlign rightClipAlign = clipFromRight(numMotifsInAllele, alignment);

    if (stutterFreeAlign.score() > leftClipAlign.score() && stutterFreeAlign.score() > rightClipAlign.score())
    {
        return stutterFreeAlign;
    }

    return (leftClipAlign.score() > rightClipAlign.score() ? leftClipAlign : rightClipAlign);
}

}
