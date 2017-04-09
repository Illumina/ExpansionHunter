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

#include <string>
using std::string;
#include <iostream>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "include/read_alignment.h"

/*****************************************************************************/

bool GetAlignFromHtsAlign(bam1_t* hts_align_ptr, Align& align,
                          bool assumeUnaligned) {
  align.name = bam_get_qname(hts_align_ptr);
  align.flag = hts_align_ptr->core.flag;

  align.status = kNoCheck;
  if (assumeUnaligned) {
    align.chrom_id = -1;       // since unaligned
    align.pos = -1;            // since unaligned
    align.mapq = 0;            // since unaligned
    align.mate_chrom_id = -1;  // since unaligned
    align.mate_pos = -1;       // since unaligned
  } else {
    align.chrom_id = hts_align_ptr->core.tid;
    align.pos = hts_align_ptr->core.pos;
    align.mapq = hts_align_ptr->core.qual;

    align.mate_chrom_id = hts_align_ptr->core.mtid;
    align.mate_pos = hts_align_ptr->core.mpos;
  }

  GetBasesFromHtsAlign(hts_align_ptr, align.bases);
  align.len = align.bases.length();
  GetQualsFromHtsAlign(hts_align_ptr, align.quals);

  return true;
}

/*****************************************************************************/

bool GetQualsFromHtsAlign(bam1_t* hts_align_ptr, string& quals) {
  uint8_t* hts_quals_ptr = bam_get_qual(hts_align_ptr);
  const int32_t read_len = hts_align_ptr->core.l_qseq;
  quals.resize(read_len);

  uint8_t* test_hts_quals_ptr = hts_quals_ptr;

  for (int32_t i = 0; i < read_len; ++i) {
    quals[i] = 33 + test_hts_quals_ptr[i];
  }

  return true;
}

/*****************************************************************************/

bool GetBasesFromHtsAlign(bam1_t* hts_align_ptr, string& bases) {
  uint8_t* hts_seq_ptr = bam_get_seq(hts_align_ptr);
  const int32_t read_len = hts_align_ptr->core.l_qseq;
  bases.resize(read_len);

  for (int32_t i = 0; i < read_len; ++i) {
    bases[i] = seq_nt16_str[bam_seqi(hts_seq_ptr, i)];
  }

  return true;
}
