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

#include "reads/Read.hh"

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
