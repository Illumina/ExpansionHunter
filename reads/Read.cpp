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

#include "reads/Read.hh"

#include <stdexcept>

using std::string;

namespace ehunter
{

Read::Read(ReadId readId, string sequence, bool isReversed)
    : readId_(std::move(readId))
    , sequence_(std::move(sequence))
    , isReversed_(isReversed)
{
    if (sequence_.empty())
    {
        std::ostringstream encoding;
        encoding << readId_;
        throw std::logic_error("Encountered empty query for " + encoding.str());
    }
}

void Read::reverseComplement()
{
    sequence_ = graphtools::reverseComplement(sequence_);
    isReversed_ = !isReversed_;
}

MappedRead::MappedRead(
    ReadId readId, string sequence, bool isReversed, int contigIndex, int64_t pos, int mapq, int mateContigIndex,
    int64_t matePos, bool isPaired, bool isMapped, bool isMateMapped)
    : Read(std::move(readId), std::move(sequence), isReversed)
    , contigIndex_(contigIndex)
    , pos_(pos)
    , mapq_(mapq)
    , mateContigIndex_(mateContigIndex)
    , matePos_(matePos)
    , isPaired_(isPaired)
    , isMapped_(isMapped)
    , isMateMapped_(isMateMapped)
{
}

bool operator==(const Read& read, const Read& mate)
{
    const bool idsAreEqual = read.readId() == mate.readId();
    const bool sequencesAreEqual = read.sequence() == mate.sequence();
    return (idsAreEqual && sequencesAreEqual);
}

bool operator==(const MappedRead& read, const MappedRead& mate)
{
    const bool readEqual = static_cast<const Read&>(read) == static_cast<const Read&>(mate);
    const bool contigEqual = read.contigIndex() == mate.contigIndex();
    const bool positionsEqual = read.pos() == mate.pos();
    const bool mapqsEqual = read.mapq() == mate.mapq();
    const bool mateContigsEqual = read.mateContigIndex() == mate.mateContigIndex();
    const bool matePositionsEqual = read.matePos() == mate.matePos();
    const bool mappingStatusesEqual = read.isMapped() == mate.isMapped();
    const bool mateMappingStatusesEqual = read.isMateMapped() == mate.isMateMapped();

    return (
        readEqual && contigEqual && positionsEqual && mapqsEqual && mateContigsEqual && matePositionsEqual
        && mappingStatusesEqual && mateMappingStatusesEqual);
}

std::ostream& operator<<(std::ostream& out, const ReadId& readId)
{
    out << readId.fragmentId() << "/" << static_cast<int>(readId.mateNumber());
    return out;
}

std::ostream& operator<<(std::ostream& out, const Read& read)
{
    out << read.readId() << " " << read.sequence();
    return out;
}

}
