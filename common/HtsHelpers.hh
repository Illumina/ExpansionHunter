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

    bool isPrimaryAlignment(bam1_t* htsAlignPtr);
    MappedRead decodeRead(bam1_t* htsAlignPtr);
    ReferenceContigInfo decodeContigInfo(bam_hdr_t* htsHeaderPtr);

} // namespace htshelpers

}
