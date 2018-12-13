//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Roman Petrovski <RPetrovski@illumina.com>
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include <algorithm>
#include <iostream>

#include "Details.hh"

namespace graphalign
{

namespace dagAligner
{

    // the 2-d table of scores filled during the alignment
    template <typename PenaltyMatrixT, bool penalizeMove> class AffineAlignMatrix
    {
    public:
        typedef PenaltyMatrixT PenaltyMatrix;

    private:
        const PenaltyMatrix penaltyMatrix_;
        const Score gapOpen_;
        const Score gapExt_;

        AlignMatrix v_;
        AlignMatrix g_;
        AlignMatrix f_;
        AlignMatrix e_;

        std::vector<typename PenaltyMatrix::QueryChar> query_;
        std::vector<typename PenaltyMatrix::TargetChar> target_;

    public:
        AffineAlignMatrix(const PenaltyMatrix& penaltyMatrix, Score gapOpen, Score gapExt)
            : penaltyMatrix_(penaltyMatrix)
            , gapOpen_(gapOpen)
            , gapExt_(gapExt)
        {
        }

        template <typename QueryIt, typename TargetIt>
        void
        init(QueryIt queryBegin, QueryIt queryEnd, TargetIt targetBegin, TargetIt targetEnd, const EdgeMap& edgeMap)
        {
            if (queryEnd == queryBegin)
            {
                throw std::logic_error("Empty query is not allowed.");
            }

            if (targetEnd == targetBegin)
            {
                throw std::logic_error("Empty target is not allowed.");
            }

            query_.clear();
            penaltyMatrix_.translateQuery(queryBegin, queryEnd, std::back_inserter(query_));
            target_.clear();
            penaltyMatrix_.translateTarget(targetBegin, targetEnd, std::back_inserter(target_));

            reset(edgeMap);

            fill(edgeMap);
        }

        typedef AlignMatrix::const_iterator const_iterator;
        template <bool localAlign> const_iterator nextBestAlign(const_iterator start, Score& bestScore) const
        {
            return !localAlign ? v_.nextBestAlign(start, queryLen() - 1, bestScore)
                               : v_.nextBestAlign(start, bestScore);
        }

        const_iterator alignBegin() const { return v_.cellOneOne(); }
        const_iterator alignEnd() const { return v_.end(); }
        int targetOffset(const_iterator cell) const { return std::distance(v_.cellOneOne(), cell) / v_.paddedRowLen(); }
        int queryOffset(const_iterator cell) const { return std::distance(v_.cellOneOne(), cell) % v_.paddedRowLen(); }
        int queryLen() const { return query_.size(); }

        bool isInsertion(int q, int t) const
        {
            const Score insExtScore = v_.at(q, t) - f_.at(q - 1, t);
            const Score insOpenScore = v_.at(q, t) - v_.at(q - 1, t);
            return gapExt_ == insExtScore || gapOpen_ + gapExt_ == insOpenScore;
        }

        bool isDeletion(int q, int t, int p) const
        {
            const Score delExtScore = v_.at(q, t) - e_.at(q, p);
            const Score delOpenScore = v_.at(q, t) - v_.at(q, p);
            return gapExt_ == delExtScore || gapOpen_ + gapExt_ == delOpenScore;
        }

        bool isMatch(int q, int t, int p) const
        {
            typename PenaltyMatrix::QueryChar queryChar = query_[q];
            typename PenaltyMatrix::TargetChar targetChar = target_[t];
            const Score alnScore = v_.at(q, t) - v_.at(q - 1, p);
            return penaltyMatrix_.isMatch(queryChar, targetChar) && penaltyMatrix_(queryChar, targetChar) == alnScore;
        }

        bool isMismatch(int q, int t, int p) const
        {
            typename PenaltyMatrix::QueryChar queryChar = query_[q];
            typename PenaltyMatrix::TargetChar targetChar = target_[t];
            const Score alnScore = v_.at(q, t) - v_.at(q - 1, p);
            return !penaltyMatrix_.isMatch(queryChar, targetChar) && penaltyMatrix_(queryChar, targetChar) == alnScore;
        }

    private:
        void computeAlignPenalties(int q, typename PenaltyMatrix::TargetChar tc, Score penalties[16])
        {
            const typename PenaltyMatrix::QueryChar* query = &query_[q];
            for (int i = 0; i < 16; ++i)
            {
                penalties[i] = penaltyMatrix_(query[i], tc);
            }
        }

        void reset(const EdgeMap& edgeMap)
        {
            const int qLen = query_.size();
            const int tLen = target_.size();

            v_.reset(qLen, tLen);
            g_.reset(qLen, tLen);
            f_.reset(qLen, tLen);
            e_.reset(qLen, tLen);

            // top left must be 0 and never change
            if (v_.at(-1, -1))
            {
                throw std::logic_error("Incorrectly initialized v_");
            }
            if (g_.at(-1, -1))
            {
                throw std::logic_error("Incorrectly initialized g_");
            }
            if (f_.at(-1, -1))
            {
                throw std::logic_error("Incorrectly initialized f_");
            }
            if (e_.at(-1, -1))
            {
                throw std::logic_error("Incorrectly initialized e_");
            }

            // first column penalises for deletion
            for (int t = 0; t < tLen; ++t)
            {
                if (penalizeMove)
                {
                    for (EdgeMap::OffsetEdges::const_iterator prevNodeIndexIt = edgeMap.prevNodesBegin(t);
                         prevNodeIndexIt != edgeMap.prevNodesEnd(t); ++prevNodeIndexIt)
                    {
                        const int p = *prevNodeIndexIt;
                        v_.at(-1, t) = std::max(v_.at(-1, t), Score(v_.at(-1, p) + gapOpen_ + gapExt_));
                        f_.at(-1, t) = std::max(v_.at(-1, t), Score(f_.at(-1, p) + gapOpen_ + gapExt_));
                    }
                }
                else
                {
                    v_.at(-1, t) = 0;
                    f_.at(-1, t) = 0;
                }
            }

            // first row penalises for insertion
            for (int q = 0; q < qLen; ++q)
            {
                v_.at(q, -1) = v_.at(q - 1, -1) + gapOpen_ + gapExt_;
                e_.at(q, -1) = e_.at(q - 1, -1) + gapOpen_ + gapExt_;
            }
        }

        void fill(const EdgeMap& edgeMap)
        {
            const int qLen = query_.size();
            const int tLen = target_.size();

            for (int t = 0; t < tLen; ++t)
            {
                for (EdgeMap::OffsetEdges::const_iterator prevNodeIndexIt = edgeMap.prevNodesBegin(t);
                     edgeMap.prevNodesEnd(t) != prevNodeIndexIt; ++prevNodeIndexIt)
                {
                    const int p = *prevNodeIndexIt;
                    recomputeForDeletion(qLen, t, p);
                    recomputeForAlign(qLen, t, p);
                }

                consolidate(qLen, t);
                recomputeForInsertion(qLen, t);
            }
        }

        void recomputeForDeletion(int qLen, int t, int p)
        {
            Score* ep = e_.row(0, p);
            Score* v = v_.row(0, p);
            Score* et = e_.row(0, t);
            for (int i = 0; i < qLen; ++i)
            {
                et[i] = std::max<Score>(et[i], std::max(ep[i] + gapExt_, v[i] + gapOpen_ + gapExt_));
            }
        }

        void recomputeForAlign(int qLen, int t, int p)
        {
            const typename PenaltyMatrix::QueryChar* query = &query_[0];

            Score* g = g_.row(0, t);
            Score* v = v_.row(-1, p);
            for (int i = 0; i < qLen; ++i)
            {
                g[i] = std::max<Score>(g[i], v[i] + penaltyMatrix_(query[i], target_[t]));
            }
        }

        void consolidate(int qLen, int t)
        {
            Score* e = e_.row(0, t);
            Score* g = g_.row(0, t);
            Score* v = v_.row(0, t);
            for (int i = 0; i < qLen; ++i)
            {
                v[i] = std::max(v[i], std::max(g[i], e[i]));
            }
        }

        void recomputeForInsertion(int qLen, int t)
        {
            Score* v = v_.row(0, t);
            Score* f = f_.row(0, t);
            Score* fp = f_.row(-1, t);
            Score* vp = v_.row(-1, t);

            for (int i = 0; i < qLen; ++i)
            {
                f[i] = std::max<Score>(f[i], std::max<Score>((fp[i] + gapExt_), (vp[i] + gapOpen_ + gapExt_)));
                v[i] = std::max(v[i], f[i]);
            }
        }

        friend std::ostream& operator<<(std::ostream& os, const AffineAlignMatrix& matrix)
        {
            return os << "AffineAlignMatrix(" << matrix.v_ << ")";
        }
    };

} // namespace dagAligner

} // namespace graphalign
