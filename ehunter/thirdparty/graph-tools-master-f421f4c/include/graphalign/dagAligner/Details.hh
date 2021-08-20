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

#include <array>
#include <cstdint>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <boost/assert.hpp>

// avoid boost lambda issues with boost::adaptors::transformed on boost 1.53
#define BOOST_RESULT_OF_USE_DECLTYPE

#include <boost/range/adaptor/transformed.hpp>

namespace graphalign
{

namespace dagAligner
{

    typedef signed short Score;
    static const Score SCORE_MIN = std::numeric_limits<Score>::min();

    /**
     * \brief Contains information about graph edges between the target sequence characters
     */
    class EdgeMap
    {
    public:
        typedef int NodeId;
        typedef std::vector<NodeId> OffsetNodeIds;
        typedef std::vector<int> OffsetEdges;

    private:
        // for target character, an id of the node to which the character belongs
        OffsetNodeIds offsetNodeIds_;

        // [index_[t], index_[t+1]) contains range in prevOffsets_ that point to base preceding t in the target graph
        std::vector<std::size_t> index_;

        OffsetEdges prevOffsets_;

    public:
        EdgeMap(
            const OffsetNodeIds& offsetNodeIds, const std::vector<std::size_t>& index,
            const std::vector<int>& prevOffsets)
            : offsetNodeIds_(offsetNodeIds)
            , index_(index)
            , prevOffsets_(prevOffsets)
        {
            if (index.back() != prevOffsets_.size())
            {
                throw std::out_of_range("index.back() != edges.size()");
            }
        }

        /**
         * \brief construct EdgeMap from list of unique node identifiers and edges as pairs of offsets
         * \param edges     offset pairs in form of {from,to} describing the connectivity in the graph
         * \param nodeIds   unique identifier of nodes in the same order as offsets appear in the edges
         * Edges have to be sorted by 'to' position and cannot create cycles ('from' < 'to')
         * Edges starting from offset -1 are ways to enter the graph. An edge from {-1,0} is implied,
         * i.e. alignments can always start at position 0.
         * For a (graph) sequence of length n, a dummy sequence of form {n,n} must be present as a marker
         * of sequence length.
         */
        EdgeMap(const std::vector<std::pair<int, int>>& edges, const std::vector<NodeId>& nodeIds)
        {
            typedef std::vector<std::pair<int, int>> Edges;
            if (edges.back().second != edges.back().first)
            {
                throw std::invalid_argument(
                    "last pair of offsets must point to itself and last character in the graph");
            }
            std::vector<NodeId>::const_iterator nodeIdIt = nodeIds.begin();

            index_.push_back(0);
            // root node offset is -1
            prevOffsets_.push_back(-1);
            index_.push_back(prevOffsets_.size());

            Edges::const_iterator edge = edges.begin();
            if (!edges.front().second)
            {
                if (-1 == edges.front().first)
                {
                    // first node edge from -1 is optional
                    ++edge;
                }
            }

            for (; edges.end() != edge;)
            {
                // fill regular offsets to previous character within the node
                while (int(offsetNodeIds_.size()) < edge->second - 1)
                {
                    prevOffsets_.push_back(int(offsetNodeIds_.size()));
                    index_.push_back(prevOffsets_.size());
                    offsetNodeIds_.push_back(*nodeIdIt);
                }

                offsetNodeIds_.push_back(*nodeIdIt);

                // insert graph edges, all incoming edges for the same target offset
                const int lastEdge = edge->second;
                do
                {
                    prevOffsets_.push_back(edge->first);
                    ++edge;
                } while (edges.end() != edge && lastEdge == edge->second);

                index_.push_back(prevOffsets_.size());

                // Now we're on a new node
                ++nodeIdIt;
            }
            // remove closing offset reference to itself
            prevOffsets_.pop_back();
            index_.pop_back();
        }

        /*
         * \param t offset of character in target character sequence
         * \return  iterator to the first offset of the predecessor character. Could be on the same or a different node
         */
        OffsetEdges::const_iterator prevNodesBegin(std::size_t t) const
        {
            assert(t < index_.size());
            return prevOffsets_.begin() + index_[t];
        }

        /*
         * \param t offset of character in target character sequence
         * \return  iterator to the one after last offset of the predecessor character. Could be on the same or a
         * different node
         */
        OffsetEdges::const_iterator prevNodesEnd(std::size_t t) const
        {
            return prevOffsets_.begin() + index_.at(t + 1);
        }

        std::size_t getNodeId(std::size_t t) const { return offsetNodeIds_.at(t); }

        friend std::ostream& operator<<(std::ostream& os, const EdgeMap& edgeMap)
        {
            return os << "EdgeMap("
                      << "nodeIds("
                      << boost::algorithm::join(
                             edgeMap.offsetNodeIds_
                                 | boost::adaptors::transformed([](std::size_t d) { return std::to_string(d); }),
                             ",")
                      << "),"
                      << "index("
                      << boost::algorithm::join(
                             edgeMap.index_
                                 | boost::adaptors::transformed([](std::size_t d) { return std::to_string(d); }),
                             ",")
                      << "),"
                      << "prevOffsets("
                      << boost::algorithm::join(
                             edgeMap.prevOffsets_
                                 | boost::adaptors::transformed([](int d) { return std::to_string(d); }),
                             ",")
                      << "))";
        }
    };

    // the 2-d table of scores filled during the alignment
    template <int pad> class PaddedAlignMatrix
    {
        std::size_t rowLen_;
        typedef std::vector<Score> Matrix;
        Matrix matrix_;

    public:
        PaddedAlignMatrix()
            : rowLen_(1)
            , matrix_(rowLen_, 0)
        {
        }

        void reset(std::size_t qLen, std::size_t tLen)
        {
            // + 1 for gap row and gap column
            rowLen_ = qLen + 1;
            matrix_.resize(paddedRowLen() * (tLen + 1));
            std::fill(matrix_.begin() + 1, matrix_.end(), SCORE_MIN);
        }

        int paddedRowLen() const { return 1 + (rowLen_ - 1 + pad - 1) / pad * pad; }
        int rowLen() const { return rowLen_; }
        Score at(int q, int t) const { return *(cellOneOne() + t * paddedRowLen() + q); }
        Score& at(int q, int t) { return *(cellOneOne() + t * paddedRowLen() + q); }
        Score* row(int q, int t) { return &(*(cellOneOne() + t * paddedRowLen() + q)); }

        typedef Matrix::iterator iterator;
        iterator cellZeroZero() { return matrix_.begin(); }
        iterator cellOneOne() { return matrix_.begin() + 1 + paddedRowLen(); }
        iterator cellOneZero() { return matrix_.begin() + 1; }
        iterator cellZeroOne() { return matrix_.begin() + paddedRowLen(); }
        typedef Matrix::const_iterator const_iterator;
        const_iterator cellOneOne() const { return matrix_.begin() + 1 + paddedRowLen(); }
        const_iterator begin() const { return matrix_.begin(); }
        const_iterator end() const { return matrix_.end() - paddedRowLen() + rowLen(); }
        const_iterator last(int q) const { return matrix_.end() - paddedRowLen() + q + 1; }
        // find best alignment after the line indicated by start iterator
        const_iterator nextBestAlign(const_iterator start, int q, Score& bestScore) const
        {
            const_iterator rowStart = begin()
                + (std::distance(begin(), start) + paddedRowLen() - 1 - q - 1) / paddedRowLen() * paddedRowLen();
            if (matrix_.end() == rowStart)
            {
                // current line is the end, no more good scores
                return end();
            }
            const_iterator rowEnd = rowStart + q + 1;
            bestScore = *rowEnd;
            const_iterator ret = rowEnd;

            while (last(q) != rowEnd)
            {
                rowEnd += paddedRowLen();
                if (bestScore < *rowEnd)
                {
                    bestScore = *rowEnd;
                    ret = rowEnd;
                }
            }
            return ret;
        }

        const_iterator nextBestAlign(const_iterator start, Score& bestScore) const
        {
            const_iterator it = start;
            if (end() == it)
            {
                // next line is the end, no more good scores
                return end();
            }
            if (!(std::distance(begin(), it) % paddedRowLen()))
            {
                ++it;
            }

            if (rowLen() == (std::distance(begin(), it) % paddedRowLen()))
            {
                it += paddedRowLen() - rowLen() + 1;
            }
            const_iterator ret = it;
            bestScore = *(it++);

            //
            while (end() != it)
            {
                if (rowLen() == (std::distance(begin(), it) % paddedRowLen()))
                {
                    it += paddedRowLen() - rowLen() + 1;
                }
                if (bestScore < *it)
                {
                    bestScore = *it;
                    ret = it;
                }
                ++it;
            }
            return ret;
        }

        int targetOffset(const_iterator cell) const { return std::distance(cellOneOne(), cell) / paddedRowLen(); }
        int queryOffset(const_iterator cell) const { return std::distance(cellOneOne(), cell) % paddedRowLen(); }

        friend std::ostream& operator<<(std::ostream& os, const PaddedAlignMatrix& matrix)
        {
            os << "AlignMatrix(\n";
            for (const_iterator it = matrix.begin(); matrix.end() > it; it += matrix.paddedRowLen())
            {
                for (std::size_t i = 0; i < matrix.rowLen_; ++i)
                {
                    os << (i ? "\t" : "[") << int(*(it + i));
                }
                os << "]\n";
            }
            return os << ")";
        }
    };

    typedef PaddedAlignMatrix<1> AlignMatrix;

    class Cigar
    {
    public:
        enum OpCode
        {
            ALIGN = 0, // 'M'
            INSERT = 1, // 'I'
            DELETE = 2, // 'D'
            SKIP = 3, // 'N' Essentially same as 'D' but not treated as a deletion.
            // Can be used for intron when aligning RNA sample against whole genome reference
            SOFT_CLIP = 4, // 'S'
            HARD_CLIP = 5, // 'H'
            PAD = 6, // 'P'
            MATCH = 7, // '='
            MISMATCH = 8, // 'X'
            UNKNOWN, // '?'
            NODE_START, // Non-standard. Indicates change of the node in graph cigar
            NODE_END // Non-standard. Indicates change of the node in graph cigar
        };

        struct Operation
        {
            typedef std::size_t ValueType;
            Operation(OpCode code, ValueType value)
                : code_(code)
                , value_(value)
            {
            }
            OpCode code_;
            ValueType value_; // normally the operation length, but for NODE contains nodeId

            char getCharCode(const char matchChar = '=', const char mismatchChar = 'X') const
            {
                const std::vector<char> CIGAR_CHARS
                    = { 'M', 'I', 'D', 'N', 'S', 'H', 'P', matchChar, mismatchChar, '?', '[', ']' };
                if (CIGAR_CHARS.size() > code_)
                {
                    return CIGAR_CHARS[code_];
                }
                throw std::out_of_range(std::string("invalid code '") + std::to_string(code_) + "'");
            }

            bool operator<(const Operation& that) const
            {
                return code_ < that.code_ || (code_ == that.code_ && value_ < that.value_);
            }

            bool operator==(const Operation& that) const { return code_ == that.code_ && value_ == that.value_; }

            friend std::ostream& operator<<(std::ostream& os, const Operation& op)
            {
                return os << "Operation(" << op.getCharCode() << op.value_ << ")";
            }
        };

        typedef std::vector<Operation> Operations;

        typedef Operations::iterator iterator;
        iterator begin() { return cigar_.begin(); }
        iterator end() { return cigar_.end(); }

        typedef Operations::const_iterator const_iterator;
        const_iterator begin() const { return cigar_.begin(); }
        const_iterator end() const { return cigar_.end(); }
        void push_back(const Operation& op) { cigar_.push_back(op); }
        Operation pop_back()
        {
            const Operation ret = cigar_.back();
            cigar_.pop_back();
            return ret;
        }
        iterator erase(const_iterator from, const_iterator to) { return cigar_.erase(from, to); }

        template <typename IteratorT> void insert(const_iterator before, IteratorT b, IteratorT e)
        {
            cigar_.insert(before, b, e);
        }

        Operation back() const { return cigar_.back(); }

        std::size_t firstNode() const
        {
            if (Cigar::NODE_START != cigar_.front().code_)
            {
                throw std::logic_error("First CIGAR op is expected to be a node start");
            }
            return cigar_.front().value_;
        }

        bool empty() const { return cigar_.empty(); }
        std::size_t length() const { return cigar_.size(); }

        OpCode lastOp() const { return cigar_.back().code_; }
        OpCode& lastOp() { return cigar_.back().code_; }
        std::size_t& lastValue() { return cigar_.back().value_; }
        void append(OpCode op)
        {
            if (lastOp() != op)
            {
                push_back(Cigar::Operation(op, 1));
            }
            else
            {
                ++lastValue();
            }
        }

        Cigar operator+(OpCode op) const
        {
            Cigar ret = *this;
            ret.append(op);
            return ret;
        }

        void reverse()
        {
            std::reverse(cigar_.begin(), cigar_.end());
            for (Operation& op : cigar_)
            {
                if (NODE_START == op.code_)
                {
                    op.code_ = NODE_END;
                }
                else if (NODE_END == op.code_)
                {
                    op.code_ = NODE_START;
                }
            }
        }

        void collapseLastEmptyNode()
        {
            if (4 < cigar_.size() && cigar_[cigar_.size() - 1].code_ == Cigar::NODE_END
                && cigar_[cigar_.size() - 4].code_ == Cigar::NODE_START
                && Cigar::SOFT_CLIP == cigar_[cigar_.size() - 3].code_
                && Cigar::DELETE == cigar_[cigar_.size() - 2].code_)
            {
                Cigar::Operation op = cigar_[cigar_.size() - 3];
                cigar_.pop_back();
                cigar_.pop_back();
                cigar_.pop_back();
                cigar_.pop_back();
                std::swap(cigar_.back(), op);
                cigar_.push_back(op);
            }
            else if (
                3 < cigar_.size() && cigar_[cigar_.size() - 1].code_ == Cigar::NODE_END
                && cigar_[cigar_.size() - 3].code_ == Cigar::NODE_START
                && Cigar::SOFT_CLIP == cigar_[cigar_.size() - 2].code_)
            {
                Cigar::Operation op = cigar_[cigar_.size() - 2];
                cigar_.pop_back();
                cigar_.pop_back();
                cigar_.pop_back();
                std::swap(cigar_.back(), op);
                cigar_.push_back(op);
            }
        }

        friend std::ostream& operator<<(std::ostream& os, const Cigar& cigar)
        {
            for (const Operation& op : cigar.cigar_)
            {
                if (NODE_START == op.code_)
                {
                    os << op.value_ << '[';
                }
                else if (NODE_END == op.code_)
                {
                    os << ']';
                }
                else
                {
                    os << op.value_ << op.getCharCode();
                }
            }
            return os;
        }

        bool operator<(const Cigar& that) const
        {
            return std::lexicographical_compare(begin(), end(), that.begin(), that.end());
        }

        bool operator==(const Cigar& that) const { return cigar_ == that.cigar_; }

    private:
        Operations cigar_;
    };

} // namespace dagAligner

} // namespace graphalign
