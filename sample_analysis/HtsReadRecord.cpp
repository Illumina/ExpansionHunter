//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "sample_analysis/HtsReadRecord.hh"

#include "common/HtsHelpers.hh"

namespace ehunter
{

HtsReadRecord::HtsReadRecord(bam1_t* htsAlignment)
    : htsAlignment_(htsAlignment)
{
}

MappedRead HtsReadRecord::decode() { return htshelpers::decodeRead(htsAlignment_); }

int HtsReadRecord::contigId() const { return htsAlignment_->core.tid; }

std::int64_t HtsReadRecord::position() const { return htsAlignment_->core.pos; }

int HtsReadRecord::mateContigId() const { return htsAlignment_->core.mtid; }

std::int64_t HtsReadRecord::matePosition() const { return htsAlignment_->core.mpos; }

}
