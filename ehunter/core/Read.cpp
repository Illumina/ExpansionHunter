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

#include "core/Read.hh"

#include <stdexcept>

using std::string;

namespace ehunter
{

bool operator==(const Read& read, const Read& mate)
{
    const bool idsAreEqual = read.readId() == mate.readId();
    const bool sequencesAreEqual = read.sequence() == mate.sequence();
    return (idsAreEqual && sequencesAreEqual);
}

bool operator==(const LinearAlignmentStats& statsA, const LinearAlignmentStats& statsB)
{
    const bool contigsEqual = statsA.chromId == statsB.chromId;
    const bool positionsEqual = statsA.pos == statsB.pos;
    const bool mapqsEqual = statsA.mapq == statsB.mapq;
    const bool mateContigsEqual = statsA.mateChromId == statsB.mateChromId;
    const bool matePositionsEqual = statsA.matePos == statsB.matePos;
    const bool mappingStatusesEqual = statsA.isMapped == statsB.isMapped;
    const bool mateMappingStatusesEqual = statsA.isMateMapped == statsB.isMateMapped;

    return (
        contigsEqual && positionsEqual && mapqsEqual && mateContigsEqual && matePositionsEqual && mappingStatusesEqual
        && mateMappingStatusesEqual);
}

std::ostream& operator<<(std::ostream& out, const LinearAlignmentStats& alignmentStats)
{
    out << "tid: " << alignmentStats.chromId << " pos: " << alignmentStats.pos
        << " mtid: " << alignmentStats.mateChromId << " mpos: " << alignmentStats.matePos
        << " mapq: " << alignmentStats.mapq << "\n";
    out << "Paired/Mapped/MateMapped: " << alignmentStats.isPaired << "/" << alignmentStats.isMapped << "/"
        << alignmentStats.isMateMapped << "\n";
    return out;
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
