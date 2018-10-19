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

#include "sample_analysis/HtsHelpers.hh"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <memory>
#include <string>

#include "thirdparty/spdlog/spdlog.h"

#include "common/sequence_operations.h"

using std::string;

namespace spd = spdlog;

namespace htshelpers
{

void DecodeQuals(bam1_t* hts_align_ptr, string& quals)
{
    uint8_t* hts_quals_ptr = bam_get_qual(hts_align_ptr);
    const int32_t readLen = hts_align_ptr->core.l_qseq;
    quals.resize(readLen);

    uint8_t* test_hts_quals_ptr = hts_quals_ptr;

    for (int32_t index = 0; index < readLen; ++index)
    {
        quals[index] = static_cast<char>(33 + test_hts_quals_ptr[index]);
    }
}

void DecodeBases(bam1_t* hts_align_ptr, string& bases)
{
    uint8_t* hts_seq_ptr = bam_get_seq(hts_align_ptr);
    const int32_t readLen = hts_align_ptr->core.l_qseq;
    bases.resize(readLen);

    for (int32_t index = 0; index < readLen; ++index)
    {
        bases[index] = seq_nt16_str[bam_seqi(hts_seq_ptr, index)];
    }
}

void DecodeAlignedRead(bam1_t* hts_align_ptr, reads::Read& read, reads::LinearAlignmentStats& alignment_stats)
{
    alignment_stats.chrom_id = hts_align_ptr->core.tid;
    alignment_stats.pos = hts_align_ptr->core.pos;
    alignment_stats.mapq = hts_align_ptr->core.qual;
    alignment_stats.mate_chrom_id = hts_align_ptr->core.mtid;
    alignment_stats.mate_pos = hts_align_ptr->core.mpos;

    uint32_t sam_flag = hts_align_ptr->core.flag;
    alignment_stats.is_mapped = sam_flag & SamFlags::kIsMapped;
    alignment_stats.is_mate_mapped = sam_flag & SamFlags::kIsMateMapped;

    read.is_first_mate = sam_flag & SamFlags::kIsFirstMate;

    const string fragment_id = bam_get_qname(hts_align_ptr);
    read.read_id = fragment_id + "/" + (read.is_first_mate ? "1" : "2");

    string bases;
    DecodeBases(hts_align_ptr, bases);
    string quals;
    DecodeQuals(hts_align_ptr, quals);

    read.sequence = lowercaseLowQualityBases(bases, quals);
}

void DecodeUnalignedRead(bam1_t* hts_align_ptr, reads::Read& read)
{
    const uint32_t sam_flag = hts_align_ptr->core.flag;
    read.is_first_mate = sam_flag & SamFlags::kIsFirstMate;

    const string fragment_id = bam_get_qname(hts_align_ptr);
    read.read_id = fragment_id + "/" + (read.is_first_mate ? "1" : "2");

    string bases;
    DecodeBases(hts_align_ptr, bases);

    string quals;
    DecodeQuals(hts_align_ptr, quals);

    read.sequence = lowercaseLowQualityBases(bases, quals);
}

} // namespace htshelpers
