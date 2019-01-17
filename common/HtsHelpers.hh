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

extern "C"
{
#include "htslib/hts.h"
#include "htslib/sam.h"
}

#include "common/ReferenceContigInfo.hh"
#include "reads/Read.hh"

namespace ehunter
{

namespace htshelpers
{

    enum SamFlags
    {
        kIsUnmapped = 0x4,
        kIsMateUnmapped = 0x8,
        kIsFirstMate = 0x40,
        kIsSecondMate = 0x80,
        // kSecondaryAlign = 0x100,
        // kSupplementaryAlignment = 0x800,
        kIsNotPrimaryLine = 0x900
    };

    LinearAlignmentStats decodeAlignmentStats(bam1_t* htsAlignPtr);
    Read decodeRead(bam1_t* htsAlignPtr);
    ReferenceContigInfo decodeContigInfo(bam_hdr_t* htsHeaderPtr);

} // namespace htshelpers

}
