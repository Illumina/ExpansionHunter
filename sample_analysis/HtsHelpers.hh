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

#include "reads/read.h"

namespace htshelpers
{

enum SamFlags
{
    kSupplementaryAlign = 0x800,
    kSecondaryAlign = 0x100,
    kIsMapped = 0x0004,
    kIsFirstMate = 0x0040,
    kIsMateMapped = 0x0008
};

void DecodeAlignedRead(bam1_t* hts_align_ptr, reads::Read& read, reads::LinearAlignmentStats& alignment_stats);
void DecodeUnalignedRead(bam1_t* hts_align_ptr, reads::Read& read);

} // namespace htshelpers