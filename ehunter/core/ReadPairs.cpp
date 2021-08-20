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

#include "core/ReadPairs.hh"

#include <stdexcept>

using std::string;
using std::vector;

namespace ehunter
{

void ReadPairs::Add(Read read)
{
    ReadPair& readPair = readPairs_[read.fragmentId()];
    const int originalMateCount = readPair.numMatesSet();

    if (read.isFirstMate() && readPair.firstMate == boost::none)
    {
        readPair.firstMate = std::move(read);
    }

    if (read.isSecondMate() && readPair.secondMate == boost::none)
    {
        readPair.secondMate = std::move(read);
    }

    const int mateCountAfterAdd = readPair.numMatesSet();

    numReads_ += mateCountAfterAdd - originalMateCount;
}

void ReadPairs::AddMateToExistingRead(Read mate)
{
    ReadPair& readPair = readPairs_.at(mate.fragmentId());
    if (mate.isFirstMate() && readPair.firstMate == boost::none)
    {
        readPair.firstMate = std::move(mate);
        ++numReads_;
    }
    else if (mate.isSecondMate() && readPair.secondMate == boost::none)
    {
        readPair.secondMate = std::move(mate);
        ++numReads_;
    }
    else
    {
        throw std::logic_error("Unable to find read placement");
    }
}

const ReadPair& ReadPairs::operator[](const string& fragment_id) const
{
    if (readPairs_.find(fragment_id) == readPairs_.end())
    {
        throw std::logic_error("Fragment " + fragment_id + " does not exist");
    }
    return readPairs_.at(fragment_id);
}

int32_t ReadPairs::NumCompletePairs() const
{
    int32_t numCompletePairs = 0;
    for (const auto& fragmentIdAndReads : readPairs_)
    {
        const ReadPair& reads = fragmentIdAndReads.second;
        if (reads.firstMate && reads.secondMate)
        {
            ++numCompletePairs;
        }
    }

    return numCompletePairs;
}

void ReadPairs::Clear()
{
    readPairs_.clear();
    numReads_ = 0;
}

bool operator==(const ReadPair& readPairA, const ReadPair& readPairB)
{
    const bool areFirstMatesEqual = readPairA.firstMate == readPairB.firstMate;
    const bool areSecondMatesEqual = readPairA.secondMate == readPairB.secondMate;
    return (areFirstMatesEqual && areSecondMatesEqual);
}

}
