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

#include <cassert>
#include <limits>
#include <ostream>

#include "graphalign/GraphAlignment.hh"

namespace ehunter
{

class StrAlign
{
public:
    enum class Type : uint8_t
    {
        kSpanning,
        kFlanking,
        kInRepeat,
        kOutside
    };

    static Type decodeType(const char typeEncoding)
    {
        switch (typeEncoding)
        {
        case 'F':
            return Type::kFlanking;
        case 'S':
            return Type::kSpanning;
        case 'I':
            return Type::kInRepeat;
        case 'O':
            return Type::kOutside;
        default:
            throw std::runtime_error("Encountered unknown StrAlign::Type: " + std::string(1, typeEncoding));
        }
    }

    StrAlign(Type type, int numMotifs, int score, int numIndels)
        : type_(type)
        , numIndels_(numIndels)
        , numMotifs_(numMotifs)
        , score_(score)
    {
        if ((numIndels < 0) or (numIndels > std::numeric_limits<decltype(numIndels_)>::max()))
        {
            throw std::runtime_error("numIndels out of range: " + std::to_string(numIndels));
        }
        if ((numMotifs < 0) or (numMotifs > std::numeric_limits<decltype(numMotifs_)>::max()))
        {
            throw std::runtime_error("numMotifs out of range: " + std::to_string(numMotifs));
        }
        if ((score < std::numeric_limits<decltype(score_)>::min())
            or (score > std::numeric_limits<decltype(score_)>::max()))
        {
            throw std::runtime_error("score out of range: " + std::to_string(score));
        }
    }

    StrAlign(char typeEncoding, int numMotifs, int score, int numIndels)
        : StrAlign(decodeType(typeEncoding), numMotifs, score, numIndels)
    {
    }

    bool operator==(const StrAlign& other) const
    {
        return type_ == other.type_ && numMotifs_ == other.numMotifs_ && score_ == other.score_
            && numIndels_ == other.numIndels_;
    }

    bool operator<(const StrAlign& other) const
    {
        if (type_ < other.type_)
            return true;
        if (type_ > other.type_)
            return false;
        if (score_ < other.score_)
            return true;
        if (score_ > other.score_)
            return false;
        if (numMotifs_ < other.numMotifs_)
            return true;
        if (numMotifs_ > other.numMotifs_)
            return false;
        return (numIndels_ < other.numIndels_);
    }

    Type type() const { return type_; }
    int numMotifs() const { return numMotifs_; }
    int score() const { return score_; }
    int numIndels() const { return numIndels_; }

private:
    Type type_;
    uint8_t numIndels_;
    uint16_t numMotifs_;
    int16_t score_;
};

std::ostream& operator<<(std::ostream& out, StrAlign::Type type);
std::ostream& operator<<(std::ostream& out, const StrAlign& summary);

class ConsistentAlignmentCalculator
{
public:
    explicit ConsistentAlignmentCalculator(int strNodeId)
        : strNodeId_(strNodeId)
    {
    }

    int strNodeId() const { return strNodeId_; }

    // Calculates longest consistent alignment by clipping from left (right)
    StrAlign clipFromLeft(int numMotifsInAllele, const graphtools::GraphAlignment& alignment) const;
    StrAlign clipFromRight(int numMotifsInAllele, const graphtools::GraphAlignment& alignment) const;

    // Calculates consistent alignment by removing PCR stutter
    StrAlign removeStutter(int numMotifsInAllele, const graphtools::GraphAlignment& alignment) const;

    StrAlign findConsistentAlignment(int numMotifsInAllele, const graphtools::GraphAlignment& alignment) const;

private:
    int matchScore_ = 5;
    int mismatchScore_ = -4;
    int gapOpenScore_ = -8;

    int strNodeId_;
};

}
