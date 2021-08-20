//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Roman Petrovski <RPetrovski@illumina.com>
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

#pragma once

#include <algorithm>
#include <iostream>
#include <string>

#include <boost/algorithm/string/join.hpp>
#include <boost/assert.hpp>

#include "dagAligner/AffineAlignMatrix.hh"
#include "dagAligner/AffineAlignMatrixVectorized.hh"
#include "dagAligner/PenaltyMatrix.hh"

namespace graphalign
{

namespace dagAligner
{
    /**
     * Performs global alignment of query against DAG of target nodes.
     * \param clipFront true instructs to represent insertions at the start of CIGAR as soft clips
     */
    template <typename AlignMatrix, bool clipFront = true> class Aligner
    {
        AlignMatrix alignMatrix_;

        // max number of best paths to backtrack
        const std::size_t maxRepeats_;

    public:
        Aligner(
            const typename AlignMatrix::PenaltyMatrix& penaltyMatrix, Score gapOpen, Score gapExt,
            std::size_t maxRepeats = 10)
            : alignMatrix_(penaltyMatrix, gapOpen, gapExt)
            , maxRepeats_(maxRepeats)
        {
        }

        template <typename QueryIt, typename TargetIt>
        void __attribute((noinline))
        align(QueryIt queryBegin, QueryIt queryEnd, TargetIt targetBegin, TargetIt targetEnd, const EdgeMap& edgeMap)
        {
            alignMatrix_.init(queryBegin, queryEnd, targetBegin, targetEnd, edgeMap);
        }

        struct Step
        {
            Cigar::OpCode operation_;
            int q_;
            int t_;
        };

        static bool removeDuplicateCigars(std::vector<Cigar>& cigars)
        {
            std::sort(cigars.begin(), cigars.end());
            std::vector<Cigar>::iterator uniq = std::unique(cigars.begin(), cigars.end());
            if (cigars.end() != uniq)
            {
                cigars.erase(uniq, cigars.end());
                return true;
            }

            return false;
        }
        template <bool localAlign>
        Score __attribute((noinline))
        backtrackAllPaths(const EdgeMap& edgeMap, std::vector<Cigar>& cigars, Score& secondBestScore) const
        {
            Score bestScore = SCORE_MIN;
            secondBestScore = SCORE_MIN;
            typename AlignMatrix::const_iterator bestCell
                = alignMatrix_.template nextBestAlign<localAlign>(alignMatrix_.alignBegin(), secondBestScore);
            for (; alignMatrix_.alignEnd() != bestCell && bestScore <= secondBestScore;
                 bestCell = alignMatrix_.template nextBestAlign<localAlign>(bestCell + 1, secondBestScore))
            {
                bestScore = secondBestScore;
                const int t = alignMatrix_.targetOffset(bestCell);
                const int q = alignMatrix_.queryOffset(bestCell);
                const int softClip = alignMatrix_.queryLen() - 1 - q;

                std::size_t firstNodeId = edgeMap.getNodeId(t);
                Cigar start;
                start.push_back(Cigar::Operation(Cigar::NODE_START, firstNodeId));
                if (softClip)
                {
                    start.push_back(Cigar::Operation(Cigar::SOFT_CLIP, softClip));
                }

                if (!backtrackPath<true>(edgeMap, start, firstNodeId, q, t, cigars))
                {
                    // ran out of cigars buffer
                    return bestScore;
                }
            }

            removeDuplicateCigars(cigars);
            if (1 < cigars.size())
            {
                // at least one duplicate, reset secondBestScore
                secondBestScore = bestScore;
            }
            else if (alignMatrix_.alignEnd() == bestCell)
            {
                // one candidate only, no second best. Might had some duplicates, reset secondBestScore
                secondBestScore = SCORE_MIN;
            }
            // else scondBest is set properly
            return bestScore;
        }

        template <bool localAlign>
        Cigar __attribute((noinline))
        backtrackBestPath(const EdgeMap& edgeMap, Score& bestScore, Score& secondBestScore) const
        {
            bestScore = SCORE_MIN;
            secondBestScore = SCORE_MIN;
            typename AlignMatrix::const_iterator bestCell
                = alignMatrix_.template nextBestAlign<localAlign>(alignMatrix_.alignBegin(), bestScore);
            if (alignMatrix_.alignEnd() == bestCell)
            {
                throw std::logic_error("No best path available");
            }
            std::vector<Cigar> ret;
            const int t = alignMatrix_.targetOffset(bestCell);
            const int q = alignMatrix_.queryOffset(bestCell);
            const int softClip = alignMatrix_.queryLen() - 1 - q;

            std::size_t firstNodeId = edgeMap.getNodeId(t);
            Cigar start;
            start.push_back(Cigar::Operation(Cigar::NODE_START, firstNodeId));
            if (softClip)
            {
                start.push_back(Cigar::Operation(Cigar::SOFT_CLIP, softClip));
            }

            backtrackPath<false>(edgeMap, start, firstNodeId, q, t, ret);
            alignMatrix_.template nextBestAlign<localAlign>(bestCell + 1, secondBestScore);
            return ret.front();
        }

    private:
        template <bool exploreAllPaths>
        Step __attribute((noinline)) stepBack(
            const EdgeMap& edgeMap, const Cigar& base, std::size_t lastNodeId, int q, int t,
            std::vector<Cigar>& cigars) const
        {
            Step ret = { Cigar::UNKNOWN, q, t };

            if (alignMatrix_.isInsertion(q, t))
            {
                const Cigar::OpCode code
                    = (1 == base.length() && Cigar::NODE_START == base.lastOp()) || Cigar::SOFT_CLIP == base.lastOp()
                    ? Cigar::SOFT_CLIP
                    : Cigar::INSERT;
                ret = Step{ code, q - 1, t };
            }

            EdgeMap::OffsetEdges::const_iterator prevNodeIndexIt = edgeMap.prevNodesBegin(t);
            while (prevNodeIndexIt != edgeMap.prevNodesEnd(t))
            {
                const int p = *prevNodeIndexIt;
                if (alignMatrix_.isDeletion(q, t, p))
                {
                    if (Cigar::UNKNOWN != ret.operation_ && exploreAllPaths)
                    {
                        // recurse if more than one path is possible
                        if (!backtrackPath<true>(edgeMap, base + Cigar::DELETE, lastNodeId, q, p, cigars))
                        {
                            return Step{ Cigar::UNKNOWN, q, t };
                        }
                    }
                    else
                    {
                        ret = Step{ Cigar::DELETE, q, p };
                    }
                }

                if (alignMatrix_.isMatch(q, t, p))
                {
                    if (Cigar::UNKNOWN != ret.operation_ && exploreAllPaths)
                    {
                        // recurse if more than one path is possible
                        if (!backtrackPath<true>(edgeMap, base + Cigar::MATCH, lastNodeId, q - 1, p, cigars))
                        {
                            return Step{ Cigar::UNKNOWN, q, t };
                        }
                    }
                    else
                    {
                        ret = Step{ Cigar::MATCH, q - 1, p };
                    }
                }
                else if (alignMatrix_.isMismatch(q, t, p))
                {
                    const Cigar::OpCode code = (1 == base.length() && Cigar::NODE_START == base.lastOp())
                            || Cigar::SOFT_CLIP == base.lastOp()
                        ? Cigar::SOFT_CLIP
                        : Cigar::MISMATCH;
                    if (Cigar::UNKNOWN != ret.operation_ && exploreAllPaths)
                    {
                        // recurse if more than one path is possible
                        if (!backtrackPath<true>(edgeMap, base + code, lastNodeId, q - 1, p, cigars))
                        {
                            return Step{ Cigar::UNKNOWN, q, t };
                        }
                    }
                    else
                    {
                        ret = Step{ code, q - 1, p };
                    }
                }

                ++prevNodeIndexIt;
            }

            if (Cigar::UNKNOWN == ret.operation_)
            {
                throw std::logic_error("backtracking failure: no nodes found on best path!");
            }

            // Return the one we have not recursed for. If only one path is possible, no recursion occurs
            return ret;
        }

        template <bool exploreAllPaths>
        bool backtrackPath(
            const EdgeMap& edgeMap, const Cigar& base, std::size_t lastNodeId, int q, int t,
            std::vector<Cigar>& cigars) const
        {
            Cigar ret = base;
            while (-1 != t)
            {
                const std::size_t curNodeId = edgeMap.getNodeId(t);
                if (lastNodeId != curNodeId)
                {
                    ret.push_back(Cigar::Operation(Cigar::NODE_END, lastNodeId));
                    ret.push_back(Cigar::Operation(Cigar::NODE_START, curNodeId));
                    lastNodeId = curNodeId;
                }

                Step step = stepBack<exploreAllPaths>(edgeMap, ret, lastNodeId, q, t, cigars);
                if (Cigar::UNKNOWN == step.operation_)
                {
                    // ran out of cigars buffer
                    return false;
                }

                ret.append(step.operation_);
                q = step.q_;
                t = step.t_;
            }

            if (-1 != q)
            {
                if (Cigar::DELETE == ret.lastOp())
                {
                    const Cigar::Operation del = ret.pop_back();
                    ret.push_back(Cigar::Operation(Cigar::INSERT, q + 1));
                    ret.push_back(del);
                }
                else
                {
                    ret.push_back(Cigar::Operation(Cigar::INSERT, q + 1));
                }
            }
            // TODO: according to while loop above, this cannot happen. Code below is legacy for
            // the original termination condition of while loop
            else if (-1 != t) // we're either on -1 row or -1 column
            {
                // count number of bases to the start of the last node
                const std::size_t nodeId = edgeMap.getNodeId(t + 1);
                int cur = t + 1;
                while (cur && edgeMap.getNodeId(cur - 1) == nodeId)
                {
                    --cur;
                }
                if (t + 1 != cur)
                {
                    ret.push_back(Cigar::Operation(Cigar::DELETE, t + 1 - cur));
                }
            }

            if (clipFront && Cigar::INSERT == ret.lastOp())
            {
                ret.lastOp() = Cigar::SOFT_CLIP;
            }
            ret.push_back(Cigar::Operation(Cigar::NODE_END, lastNodeId));

            ret.collapseLastEmptyNode();
            ret.reverse();
            if (cigars.size() == maxRepeats_ && !removeDuplicateCigars(cigars))
            {
                return false;
            }
            cigars.push_back(ret);
            return true;
        }

        friend std::ostream& operator<<(std::ostream& os, const Aligner& aligner)
        {
            return os << "Aligner(" << aligner.alignMatrix_ << ")";
        }
    };

} // namespace dagAligner

// template <bool penalizeMove>
// class DagAligner
//     : public dagAligner::Aligner<dagAligner::AffineAlignMatrix<dagAligner::FixedPenaltyMatrix, penalizeMove>>
// {
// public:
// DagAligner(const dagAligner::FixedPenaltyMatrix& penaltyMatrix, dagAligner::Score gapOpen, dagAligner::Score gapExt)
// : dagAligner::Aligner<dagAligner::AffineAlignMatrix<dagAligner::FixedPenaltyMatrix, penalizeMove>>(
// penaltyMatrix, gapOpen, gapExt)
// {
// }
//
// DagAligner(dagAligner::Score match, dagAligner::Score mismatch, dagAligner::Score gapOpen, dagAligner::Score gapExt)
// : dagAligner::Aligner<dagAligner::AffineAlignMatrix<dagAligner::FixedPenaltyMatrix, penalizeMove>>(
// dagAligner::FixedPenaltyMatrix(match, mismatch), gapOpen, gapExt)
// {
// }
// };

template <bool penalizeMove, bool clipFront = true, bool matchQueryN = true, bool matchTargetN = true>
class DagAligner : public dagAligner::Aligner<
                       dagAligner::AffineAlignMatrixVectorized<
                           dagAligner::FixedPenaltyMatrix<matchQueryN, matchTargetN>, penalizeMove>,
                       clipFront>
{
    typedef dagAligner::FixedPenaltyMatrix<matchQueryN, matchTargetN> PenaltyMatrix;

public:
    DagAligner(const PenaltyMatrix& penaltyMatrix, dagAligner::Score gapOpen, dagAligner::Score gapExt)
        : dagAligner::Aligner<dagAligner::AffineAlignMatrixVectorized<PenaltyMatrix, penalizeMove>, clipFront>(
              penaltyMatrix, gapOpen, gapExt)
    {
    }

    DagAligner(dagAligner::Score match, dagAligner::Score mismatch, dagAligner::Score gapOpen, dagAligner::Score gapExt)
        : dagAligner::Aligner<dagAligner::AffineAlignMatrixVectorized<PenaltyMatrix, penalizeMove>, clipFront>(
              PenaltyMatrix(match, mismatch), gapOpen, gapExt)
    {
    }
};

// template <bool penalizeMove>
// class DagAligner
//     : public dagAligner::Aligner<dagAligner::AffineAlignMatrix<dagAligner::FreePenaltyMatrix, penalizeMove>>
// {
// public:
//    DagAligner(const dagAligner::FreePenaltyMatrix& penaltyMatrix, dagAligner::Score gapOpen, dagAligner::Score
//    gapExt)
//        : dagAligner::Aligner<dagAligner::AffineAlignMatrix<dagAligner::FreePenaltyMatrix, penalizeMove>>(
//              penaltyMatrix, gapOpen, gapExt)
//    {
//    }
//
//    DagAligner(dagAligner::Score match, dagAligner::Score mismatch, dagAligner::Score gapOpen, dagAligner::Score
//    gapExt)
//        : dagAligner::Aligner<dagAligner::AffineAlignMatrix<dagAligner::FreePenaltyMatrix, penalizeMove>>(
//              dagAligner::FreePenaltyMatrix(match, mismatch), gapOpen, gapExt)
//    {
//    }
// };

} // namespace graphalign
