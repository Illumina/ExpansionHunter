//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include <map>
#include <ostream>
#include <string>
#include <vector>

namespace ehunter
{

enum class ReadType
{
    kSpanning,
    kFlanking,
    kRepeat,
    kOther
};

enum class Sex
{
    kMale,
    kFemale
};

Sex decodeSampleSex(const std::string& encoding);

enum class AlleleCount
{
    kZero,
    kOne,
    kTwo
};

class NumericInterval
{
public:
    NumericInterval()
        : start_(0)
        , end_(0)
    {
    }

    NumericInterval(int start, int end)
        : start_(start)
        , end_(end)
    {
    }

    int start() const { return start_; }
    int end() const { return end_; }

    NumericInterval& operator=(const NumericInterval& other)
    {
        start_ = other.start_;
        end_ = other.end_;
        return *this;
    }

    bool operator==(const NumericInterval& other) const { return start_ == other.start_ && end_ == other.end_; }

private:
    int start_;
    int end_;
};

template <typename T> struct LabeledSequence
{
    LabeledSequence(const std::string sequence, T label)
        : sequence(sequence)
        , label(label)
    {
    }

    bool operator==(const LabeledSequence<T>& other) const
    {
        return sequence == other.sequence && label == other.label;
    }

    std::string sequence;
    T label;
};

std::ostream& operator<<(std::ostream& out, ReadType readType);
std::ostream& operator<<(std::ostream& out, AlleleCount alleleCount);
std::ostream& operator<<(std::ostream& out, NumericInterval numericInterval);

}
