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

#include "reads/ReadPairs.hh"

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
        if (reads.firstMate != boost::none && reads.secondMate != boost::none)
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
