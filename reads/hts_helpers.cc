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
#include <memory>
#include <string>

#include "third_party/spdlog/spdlog.h"

using std::string;

namespace spd = spdlog;

namespace htshelpers {

void DecodeQuals(bam1_t* hts_align_ptr, string& quals) {
  uint8_t* hts_quals_ptr = bam_get_qual(hts_align_ptr);
  const int32_t read_len = hts_align_ptr->core.l_qseq;
  quals.resize(read_len);

  uint8_t* test_hts_quals_ptr = hts_quals_ptr;

  for (int32_t index = 0; index < read_len; ++index) {
    quals[index] = static_cast<char>(33 + test_hts_quals_ptr[index]);
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

void DecodeAlignedRead(bam1_t* hts_align_ptr, reads::ReadPtr& read_ptr) {
  const string name = bam_get_qname(hts_align_ptr);

  string bases;
  DecodeBases(hts_align_ptr, bases);

  string quals;
  DecodeQuals(hts_align_ptr, quals);

  read_ptr = std::make_shared<reads::Read>(name, bases, quals);

  read_ptr->SetSamChromId(hts_align_ptr->core.tid);
  read_ptr->SetSamPos(hts_align_ptr->core.pos);
  read_ptr->SetSamMapq(hts_align_ptr->core.qual);
  read_ptr->SetSamMateChromId(hts_align_ptr->core.mtid);
  read_ptr->SetSamMatePos(hts_align_ptr->core.mpos);

  uint32_t sam_flag = hts_align_ptr->core.flag;
  read_ptr->SetIsSamMapped(sam_flag & SamFlags::kIsMapped);
  read_ptr->SetIsFirstMate(sam_flag & SamFlags::kIsFirstMate);
  read_ptr->SetIsMateSamMapped(sam_flag & SamFlags::kIsMateMapped);
}

void DecodeUnalignedRead(bam1_t* hts_align_ptr, reads::ReadPtr& read_ptr) {
  const string name = bam_get_qname(hts_align_ptr);

  string bases;
  DecodeBases(hts_align_ptr, bases);

  string quals;
  DecodeQuals(hts_align_ptr, quals);

  read_ptr = std::make_shared<reads::Read>(name, bases, quals);

  read_ptr->SetSamChromId(-1);
  read_ptr->SetIsSamMapped(false);
  read_ptr->SetSamPos(-1);
  read_ptr->SetSamMapq(0);
  read_ptr->SetIsMateSamMapped(false);
  read_ptr->SetSamMateChromId(-1);
  read_ptr->SetSamMatePos(-1);
}

}  // namespace htshelpers