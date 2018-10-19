// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.

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
    template <typename PenaltyMatrixT, bool penalizeMove, int step = 16> class AffineAlignMatrixVectorized
    {
    public:
        typedef PenaltyMatrixT PenaltyMatrix;
    private:
        const PenaltyMatrix penaltyMatrix_;
        const Score gapOpen_;
        const Score gapExt_;

        PaddedAlignMatrix<step> v_;
        PaddedAlignMatrix<step> g_;
        PaddedAlignMatrix<step> f_;
        PaddedAlignMatrix<step> e_;

        std::vector<typename PenaltyMatrix::QueryChar> query_;
        std::vector<typename PenaltyMatrix::TargetChar> target_;
        std::vector<Score> alignmentPenalties_[PenaltyMatrix::TARGET_CHAR_MAX_ + 1];

    public:
        AffineAlignMatrixVectorized(const PenaltyMatrix& penaltyMatrix, Score gapOpen, Score gapExt)
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

            // avoid "uninitialized read" complaints from valgrind
            query_.resize((std::distance(queryBegin, queryEnd) + step - 1) / step * step);
            query_.clear();
            penaltyMatrix_.translateQuery(queryBegin, queryEnd, std::back_inserter(query_));
            target_.clear();
            penaltyMatrix_.translateTarget(targetBegin, targetEnd, std::back_inserter(target_));

            reset(edgeMap);

            fill(edgeMap);
            // std::cerr << *this << std::endl;
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
            // std::cerr << "Score " << v_.at(q, t) << std::endl;
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
            return penaltyMatrix_.isMatch(queryChar, targetChar) && penaltyMatrix_(queryChar, queryChar) == alnScore;
        }

        bool isMismatch(int q, int t, int p) const
        {
            typename PenaltyMatrix::QueryChar queryChar = query_[q];
            typename PenaltyMatrix::TargetChar targetChar = target_[t];
            const Score alnScore = v_.at(q, t) - v_.at(q - 1, p);
            return !penaltyMatrix_.isMatch(queryChar, targetChar) && penaltyMatrix_(queryChar, targetChar) == alnScore;
        }

    private:

        // __attribute((noinline))
        void reset(const EdgeMap& edgeMap)
        {
            const int qLen = query_.size();
            const int tLen = target_.size();

            for (typename PenaltyMatrix::TargetChar tc = 0; tc <= PenaltyMatrix::TARGET_CHAR_MAX_; ++tc)
            {
                alignmentPenalties_[tc].resize((qLen + step - 1) / step * step, 0);
                for (int q = 0; q < qLen; q += step)
                {
                    computeAlignPenalties(q, tc, &alignmentPenalties_[tc].front() + q);
                }
            }

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

            // first row penalizes for insertion
            for (int q = 0; q < queryLen(); ++q)
            {
                v_.at(q, -1) = v_.at(q - 1, -1) + gapOpen_ + gapExt_;
                e_.at(q, -1) = e_.at(q - 1, -1) + gapOpen_ + gapExt_;
            }
        }

        // __attribute((noinline))
        void fill(const EdgeMap& edgeMap)
        {
            const int qLen = query_.size();
            const int tLen = target_.size();

            for (int t = 0; t < tLen; ++t)
            {
                const typename PenaltyMatrix::TargetChar tc = target_[t];
                const Score* const penalties = &alignmentPenalties_[tc].front();
                for (EdgeMap::OffsetEdges::const_iterator prevNodeIndexIt = edgeMap.prevNodesBegin(t);
                     edgeMap.prevNodesEnd(t) != prevNodeIndexIt; ++prevNodeIndexIt)
                {
                    const int p = *prevNodeIndexIt;
                    for (int q = 0; q < qLen; q += step)
                    {
                        recomputeForDeletion(q, t, p);
                        recomputeForAlign(q, t, p, penalties + q);
                    }
                }

                for (int q = 0; q < qLen; q += step)
                {
                    consolidate(q, t);
                }
                recomputeForInsertion(qLen, t);
            }
        }

        // __attribute((noinline))
        void computeAlignPenalties(int q, typename PenaltyMatrix::TargetChar tc, Score penalties[step])
        {
            const typename PenaltyMatrix::QueryChar* query = &query_[q];
            for (int i = 0; i < step; ++i)
            {
                penalties[i] = penaltyMatrix_(query[i], tc);
            }
        }

        // __attribute((noinline))
        void recomputeForAlign(int q, int t, int p, const Score penalties[step])
        {
            Score tmp[step];
            Score* v = v_.row(q - 1, p);
            for (int i = 0; i < step; ++i)
            {
                tmp[i] = v[i] + penalties[i];
            }

            Score* g = g_.row(q, t);
            for (int i = 0; i < step; ++i)
            {
                g[i] = g[i] > tmp[i] ? g[i] : tmp[i];
            }
        }

        // __attribute((noinline))
        void recomputeForDeletion(int q, int t, int p)
        {
            Score* ep = e_.row(q, p);
            Score tmpEp[step];
            for (int i = 0; i < step; ++i)
            {
                tmpEp[i] = ep[i] + gapExt_;
            }

            Score* vp = v_.row(q, p);
            Score tmpVp[step];
            for (int i = 0; i < step; ++i)
            {
                tmpVp[i] = vp[i] + gapOpen_ + gapExt_;
            }

            Score tmp[step];
            for (int i = 0; i < step; ++i)
            {
                tmp[i] = tmpEp[i] > tmpVp[i] ? tmpEp[i] : tmpVp[i];
            }

            Score* e = e_.row(q, t);
            for (int i = 0; i < step; ++i)
            {
                e[i] = e[i] > tmp[i] ? e[i] : tmp[i];
            }
        }

        // __attribute((noinline))
        void recomputeForInsertion(int qLen, int t)
        {
            Score* v = v_.row(0, t);
            Score* f = f_.row(0, t);
            Score* fp = f_.row(-1, t);
            Score* vp = v_.row(-1, t);

            // cannot be vectorized since insertion is computed horisontally
            for (int i = 0; i < qLen; ++i)
            {
                f[i] = std::max<Score>(f[i], std::max<Score>((fp[i] + gapExt_), (vp[i] + gapOpen_ + gapExt_)));
                v[i] = std::max(v[i], f[i]);
            }
        }

        // __attribute((noinline))
        void consolidate(int q, int t)
        {
            Score tmp[step];
            Score* e = e_.row(q, t);
            Score* g = g_.row(q, t);
            for (int i = 0; i < step; ++i)
            {
                tmp[i] = g[i] > e[i] ? g[i] : e[i];
            }

            Score* v = v_.row(q, t);
            for (int i = 0; i < step; ++i)
            {
                v[i] = v[i] > tmp[i] ? v[i] : tmp[i];
            }
        }

        friend std::ostream& operator<<(std::ostream& os, const AffineAlignMatrixVectorized& matrix)
        {
            return os << "AffineAlignMatrix(" << matrix.v_ << ")";
        }
    };

} // namespace dagAligner

} // namespace graphalign
