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
#include <ostream>

#include "graphalign/GraphAlignment.hh"

namespace ehunter
{

class StrAlign
{
public:
    enum class Type
    {
        kSpanning,
        kFlanking,
        kInRepeat,
        kOutside
    };

    StrAlign(Type type, int numMotifs, int score, int numIndels)
        : type_(type)
        , numMotifs_(numMotifs)
        , score_(score)
        , numIndels_(numIndels)
    {
    }

    StrAlign(char typeEncoding, int numMotifs, int score, int numIndels)
        : type_(Type::kOutside)
        , numMotifs_(numMotifs)
        , score_(score)
        , numIndels_(numIndels)
    {
        switch (typeEncoding)
        {
        case 'F':
            type_ = Type::kFlanking;
            break;
        case 'S':
            type_ = Type::kSpanning;
            break;
        case 'I':
            type_ = Type::kInRepeat;
            break;
        case 'O':
            type_ = Type::kOutside;
            break;
        default:
            std::string message = "Encountered unknown StrAlign::Type: ";
            message.push_back(typeEncoding);
            throw std::runtime_error(message);
        }
    }

    bool operator==(const StrAlign& other) const
    {
        return type_ == other.type_ && numMotifs_ == other.numMotifs_ && score_ == other.score_
            && numIndels_ == other.numIndels_;
    }

    Type type() const { return type_; }
    int numMotifs() const { return numMotifs_; }
    int score() const { return score_; }
    int numIndels() const { return numIndels_; }

private:
    Type type_;
    int numMotifs_;
    int score_;
    int numIndels_;
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
    StrAlign clipFromLeft(int numMotifsInAllele, const graphtools::GraphAlignment& alignment);
    StrAlign clipFromRight(int numMotifsInAllele, const graphtools::GraphAlignment& alignment);

    // Calculates consistent alignment by removing PCR stutter
    StrAlign removeStutter(int numMotifsInAllele, const graphtools::GraphAlignment& alignment);

    StrAlign findConsistentAlignment(int numMotifsInAllele, const graphtools::GraphAlignment& alignment);

private:
    int matchScore_ = 5;
    int mismatchScore_ = -4;
    int gapOpenScore_ = -8;

    int strNodeId_;
};

}
