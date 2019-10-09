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

#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "graphalign/GraphAlignment.hh"
#include "graphalign/LinearAlignmentParameters.hh"
#include "graphcore/Graph.hh"

namespace ehunter
{

// Summarizes information pertaining to one alignment of a read to an STR.
// Note that different ways of soft-clipping the read count as different alignments.
class StrAlignment
{
public:
    enum class Type
    {
        kSpanning,
        kFlanking,
        kInrepeat
    };

    /**
     * @param numUnits: Number of repeat units overlapped by this alignment.
     * @param readType Used to distinguish ReadType::kSpanning from non-spanning reads.
     * @param score: Log-probability of observing the read if numUnits represents the true allele.
     * @param clippedReadLength: length of the aligned portion of the read
     *
     * The probability (score) need not be normalized, but the normalization constant should be the same for all scores.
     * This means that alignment scores can be used directly, provided the same scoring system is used throughout.
     * Any clipped bases should already have been penalized (e.g for normalized scores, add log(1/4) per clipped base).
     * The penalty per clipped base should be the same as randomBasePenalty, which for now (when using the alignment
     * scores directly) we are recommending to set to 0. TODO: should we use (5+log(1/4)) instead?
     */
    StrAlignment(int numUnits, Type type, int score, int clippedReadLength)
        : numUnits_(numUnits)
        , type_(type)
        , score_(score)
        , clippedReadLength_(clippedReadLength)
    {
    }

    bool operator==(const StrAlignment& other) const
    {
        return type_ == other.type_ && numUnits_ == other.numUnits_ && score_ == other.score_;
    }

    bool operator<(const StrAlignment& other) const
    {
        if (type_ != other.type_)
        {
            return type_ < other.type_;
        }

        if (numUnits_ != other.numUnits_)
        {
            return numUnits_ < other.numUnits_;
        }

        return score_ < other.score_;
    }

    int numUnits() const { return numUnits_; }
    Type type() const { return type_; }
    int score() const { return score_; }
    int clippedReadLength() const { return clippedReadLength_; }

    bool isSpanning() const { return type_ == Type::kSpanning; }
    bool isRepeat() const { return type_ == Type::kInrepeat; }

private:
    int numUnits_;
    Type type_;
    int score_;
    int clippedReadLength_;
};

std::ostream& operator<<(std::ostream& os, const StrAlignment::Type& alignmentType);
std::ostream& operator<<(std::ostream& os, const StrAlignment& alignmentSummary);

// Summarizes information pertaining to all high-scoring alignments of one read to an STR
class ReadSummaryForStr
{
public:
    explicit ReadSummaryForStr(int readLength)
        : readLength_(readLength)
    {
    }

    int readLength() const { return readLength_; }
    bool hasAlignments() const { return !alignments_.empty(); }
    int numAlignments() const { return alignments_.size(); }
    const std::vector<StrAlignment>& alignments() const { return alignments_; }
    void addAlignment(StrAlignment alignment) { alignments_.push_back(std::move(alignment)); }

    bool operator==(const ReadSummaryForStr& other) const
    {
        return readLength_ == other.readLength_ && alignments_ == other.alignments_;
    }

private:
    int readLength_;
    std::vector<StrAlignment> alignments_;
};

// Summarizes information pertaining to one alignment of a read to a small variant.
class SmallVariantAlignment
{
public:
    enum class Type
    {
        kSpanning,
        kUpstreamFlanking,
        kDownstreamFlanking
    };

    SmallVariantAlignment(graphtools::NodeId nodeId, Type type, int score)
        : nodeId_(nodeId)
        , type_(type)
        , score_(score)
    {
    }

    bool operator==(const SmallVariantAlignment& other) const
    {
        return nodeId_ == other.nodeId_ && type_ == other.type_ && score_ == other.score_;
    }

    bool operator<(const SmallVariantAlignment& other) const
    {
        if (nodeId_ != other.nodeId_)
        {
            return nodeId_ < other.nodeId_;
        }

        if (type_ != other.type_)
        {
            return type_ < other.type_;
        }

        return score_ < other.score_;
    }

    graphtools::NodeId nodeId() const { return nodeId_; }
    Type type() const { return type_; }
    int score() const { return score_; }

private:
    graphtools::NodeId nodeId_;
    Type type_;
    int score_;
};

std::ostream& operator<<(std::ostream& os, const SmallVariantAlignment::Type& alignmentType);
std::ostream& operator<<(std::ostream& os, const SmallVariantAlignment& alignmentSummary);

// Summarizes information pertaining to all high-scoring alignments of one read to a small variant
class ReadSummaryForSmallVariant
{
public:
    explicit ReadSummaryForSmallVariant(int readLength)
        : readLength_(readLength)
    {
    }

    int readLength() const { return readLength_; }
    bool hasAlignments() const { return !alignments_.empty(); }
    int numAlignments() const { return alignments_.size(); }
    const std::vector<SmallVariantAlignment>& alignments() const { return alignments_; }
    void addAlignment(SmallVariantAlignment alignment) { alignments_.push_back(alignment); }

    bool operator==(const ReadSummaryForSmallVariant& other) const
    {
        return readLength_ == other.readLength_ && alignments_ == other.alignments_;
    }

private:
    int readLength_;
    std::vector<SmallVariantAlignment> alignments_;
};

int scoreAlignment(
    const graphtools::GraphAlignment& alignment, LinearAlignmentParameters parameters = LinearAlignmentParameters());

}
