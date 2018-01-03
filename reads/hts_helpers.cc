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

#include "reads/hts_helpers.h"

#include <cstdint>
#include <string>

using std::string;

namespace htshelpers {

void DecodeQuals(bam1_t* hts_align_ptr, string& quals) {
  uint8_t* hts_quals_ptr = bam_get_qual(hts_align_ptr);
  const int32_t read_len = hts_align_ptr->core.l_qseq;
  quals.resize(read_len);

  uint8_t* test_hts_quals_ptr = hts_quals_ptr;

  for (int32_t index = 0; index < read_len; ++index) {
    quals[index] = 33 + test_hts_quals_ptr[index];
  }
}

void DecodeBases(bam1_t* hts_align_ptr, string& bases) {
  uint8_t* hts_seq_ptr = bam_get_seq(hts_align_ptr);
  const int32_t read_len = hts_align_ptr->core.l_qseq;
  bases.resize(read_len);

  for (int32_t index = 0; index < read_len; ++index) {
    bases[index] = seq_nt16_str[bam_seqi(hts_seq_ptr, index)];
  }
}

void DecodeAlignedRead(bam1_t* hts_align_ptr, reads::Read& read) {
  const string name = bam_get_qname(hts_align_ptr);

  string bases;
  DecodeBases(hts_align_ptr, bases);

  string quals;
  DecodeQuals(hts_align_ptr, bases);

  read.SetCoreInfo(name, bases, quals);

  read.SetSamChromId(hts_align_ptr->core.tid);
  read.SetSamPos(hts_align_ptr->core.pos);
  read.SetSamMapq(hts_align_ptr->core.qual);
  read.SetSamMateChromId(hts_align_ptr->core.mtid);
  read.SetSamMatePos(hts_align_ptr->core.mpos);

  uint32_t sam_flag = hts_align_ptr->core.flag;
  read.SetIsSamMapped(sam_flag & SamFlags::kIsMapped);
  read.SetIsFirstMate(sam_flag & SamFlags::kIsFirstMate);
  read.SetIsMateSamMapped(sam_flag & SamFlags::kIsMateMapped);
}

void DecodeUnalignedRead(bam1_t* hts_align_ptr, reads::Read& read) {
  const string name = bam_get_qname(hts_align_ptr);

  string bases;
  DecodeBases(hts_align_ptr, bases);

  string quals;
  DecodeQuals(hts_align_ptr, bases);

  read.SetCoreInfo(name, bases, quals);

  read.SetSamChromId(-1);
  read.SetIsSamMapped(false);
  read.SetSamPos(-1);
  read.SetSamMapq(0);
  read.SetIsMateSamMapped(false);
  read.SetSamMateChromId(-1);
  read.SetSamMatePos(-1);
}

}  // namespace htshelpers