//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#pragma once

#include <cstdint>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/functional/hash.hpp>

#include "classification/AlignmentClassifier.hh"
#include "graphalign/GraphAlignment.hh"

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
    Read(ReadId readId, std::string sequence)
        : readId_(std::move(readId))
        , sequence_(std::move(sequence))
    {
        if (sequence_.empty())
        {
            std::ostringstream encoding;
            encoding << readId_;
            throw std::logic_error("Encountered empty query for " + encoding.str());
        }
    }

    const ReadId& readId() const { return readId_; }
    const FragmentId& fragmentId() const { return readId_.fragmentId(); }
    MateNumber mateNumber() const { return readId_.mateNumber(); }
    const std::string& sequence() const { return sequence_; }
    void setSequence(std::string sequence) { sequence_ = sequence; }

    bool isFirstMate() const { return mateNumber() == MateNumber::kFirstMate; }
    bool isSecondMate() const { return mateNumber() == MateNumber::kSecondMate; }

private:
    ReadId readId_;
    std::string sequence_;
};

struct LinearAlignmentStats
{
    int32_t chromId = -1;
    int32_t pos = -1;
    int32_t mapq = -1;
    int32_t mateChromId = -1;
    int32_t matePos = -1;
    bool isMapped = false;
    bool isMateMapped = false;
};

using ReadIdToLinearAlignmentStats = std::unordered_map<std::string, LinearAlignmentStats>;

bool operator==(const Read& read, const Read& mate);
bool operator==(const LinearAlignmentStats& statsA, const LinearAlignmentStats& statsB);

class RepeatAlignmentStats
{
public:
    RepeatAlignmentStats(
        const GraphAlignment& canonical_alignment, AlignmentType canonical_alignment_type,
        int32_t num_repeat_units_spanned)
        : canonical_alignment_(canonical_alignment)
        , canonical_alignment_type_(canonical_alignment_type)
        , num_repeat_units_spanned_(num_repeat_units_spanned)
    {
    }

    const GraphAlignment& canonicalAlignment() const { return canonical_alignment_; }
    AlignmentType canonicalAlignmentType() const { return canonical_alignment_type_; }
    int32_t numRepeatUnitsSpanned() const { return num_repeat_units_spanned_; }

private:
    GraphAlignment canonical_alignment_;
    AlignmentType canonical_alignment_type_;
    int32_t num_repeat_units_spanned_;
};

using ReadIdToRepeatAlignmentStats = std::unordered_map<std::string, RepeatAlignmentStats>;

std::ostream& operator<<(std::ostream& out, const Read& read);

}
