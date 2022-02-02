//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include <map>
#include <ostream>
#include <sstream>
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

enum class ChromType
{
    kX,
    kY,
    kAutosome
};

Sex decodeSampleSex(const std::string& encoding);

enum class AlleleCount
{
    kOne = 1,
    kTwo = 2
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

template <typename T> std::string streamToString(const T& streamable)
{
    std::stringstream out;
    out << streamable;
    return out.str();
}

std::ostream& operator<<(std::ostream& out, Sex sex);
std::ostream& operator<<(std::ostream& out, ReadType readType);
std::ostream& operator<<(std::ostream& out, AlleleCount alleleCount);
std::ostream& operator<<(std::ostream& out, NumericInterval numericInterval);

/// \brief Returns true if the path refers to a URL instead of a local file
///
/// This does not test if the URL is well formed
///
bool isURL(const std::string& path);

}
