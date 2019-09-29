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

#include <cstdint>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/functional/hash.hpp>

#include "graphalign/GraphAlignment.hh"
#include "graphutils/SequenceOperations.hh"

namespace ehunter
{

using FragmentId = std::string;

enum class MateNumber
{
    kFirstMate = 1,
    kSecondMate = 2
};

class ReadId
{
public:
    ReadId(FragmentId fragmentId, MateNumber mateNumber)
        : fragmentId_(std::move(fragmentId))
        , mateNumber_(mateNumber)
    {
        if (fragmentId_.empty())
        {
            throw std::logic_error("Encountered an empty fragment id");
        }
    }

    const FragmentId& fragmentId() const { return fragmentId_; }
    MateNumber mateNumber() const { return mateNumber_; }

    bool operator==(const ReadId& other) const
    {
        return fragmentId_ == other.fragmentId_ && mateNumber_ == other.mateNumber_;
    }

    friend std::size_t hash_value(const ReadId& readId)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, readId.fragmentId_);
        boost::hash_combine(seed, static_cast<int>(readId.mateNumber_));

        return seed;
    }

private:
    FragmentId fragmentId_;
    MateNumber mateNumber_;
};

std::ostream& operator<<(std::ostream& out, const ReadId& readId);

class Read
{
public:
    Read(ReadId readId, std::string sequence, bool isReversed);
    const ReadId& readId() const { return readId_; }
    const FragmentId& fragmentId() const { return readId_.fragmentId(); }
    MateNumber mateNumber() const { return readId_.mateNumber(); }
    const std::string& sequence() const { return sequence_; }

    bool isFirstMate() const { return mateNumber() == MateNumber::kFirstMate; }
    bool isSecondMate() const { return mateNumber() == MateNumber::kSecondMate; }
    // Return whether the read is reverse complemented relative to it's
    //  original direction during sequencing
    bool isReversed() const { return isReversed_; }
    void reverseComplement();

protected:
    ReadId readId_;
    std::string sequence_;
    bool isReversed_;
};

class MappedRead : public Read
{
public:
    MappedRead(
        ReadId readId, std::string sequence, bool isReversed, int contigIndex, int64_t pos, int mapq,
        int mateContigIndex, int64_t matePos, bool isPaired, bool isMapped, bool isMateMapped);

    int contigIndex() const { return contigIndex_; }
    int64_t pos() const { return pos_; }
    int64_t approximateEnd() const { return pos() + sequence().length(); }
    int mapq() const { return mapq_; }
    int mateContigIndex() const { return mateContigIndex_; }
    int64_t matePos() const { return matePos_; }
    bool isPaired() const { return isPaired_; }
    bool isMapped() const { return isMapped_; }
    bool isMateMapped() const { return isMateMapped_; }

private:
    int contigIndex_;
    int64_t pos_;
    int mapq_;
    int mateContigIndex_;
    int64_t matePos_;
    bool isPaired_;
    bool isMapped_;
    bool isMateMapped_;
};

class ReadRecordWrapper
{
};

bool operator==(const Read& read, const Read& mate);
bool operator==(const MappedRead& read, const MappedRead& mate);

std::ostream& operator<<(std::ostream& out, const Read& read);

}
