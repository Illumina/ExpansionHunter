//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "io/SampleStats.hh"

#include <cstdint>
#include <string>
#include <utility>

#include <boost/filesystem.hpp>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "core/HtsHelpers.hh"

using std::pair;
using std::string;
using std::vector;

namespace ehunter
{

int extractReadLength(const string& htsFilePath)
{
    samFile* htsFilePtr = sam_open(htsFilePath.c_str(), "r");
    if (!htsFilePtr)
    {
        throw std::runtime_error("Failed to read " + htsFilePath);
    }
    bam_hdr_t* htsHeaderPtr = sam_hdr_read(htsFilePtr);
    if (!htsHeaderPtr)
    {
        throw std::runtime_error("Failed to read the header of " + htsFilePath);
    }

    enum
    {
        kSupplementaryAlign = 0x800,
        kSecondaryAlign = 0x100
    };

    int readLength = 99;
    bam1_t* htsAlignmentPtr = bam_init1();
    int ret;
    while ((ret = sam_read1(htsFilePtr, htsHeaderPtr, htsAlignmentPtr)) >= 0)
    {
        const bool isSupplementary = htsAlignmentPtr->core.flag & kSupplementaryAlign;
        const bool isSecondary = htsAlignmentPtr->core.flag & kSecondaryAlign;
        const bool isPrimaryAlign = (!isSupplementary) && (!isSecondary);
        if (isPrimaryAlign)
        {
            readLength = htsAlignmentPtr->core.l_qseq;
            break;
        }
    }

    if (ret < 0)
    {
        throw std::runtime_error("Failed to extract a read from " + htsFilePath);
    }

    bam_destroy1(htsAlignmentPtr);
    bam_hdr_destroy(htsHeaderPtr);
    sam_close(htsFilePtr);

    return readLength;
}

ReferenceContigInfo extractReferenceContigInfo(const std::string& htsFilePath)
{
    std::unique_ptr<samFile, decltype(&hts_close)> htsFilePtr(hts_open(htsFilePath.c_str(), "r"), hts_close);
    if (!htsFilePtr)
    {
        throw std::runtime_error("Failed to read " + htsFilePath);
    }

    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> htsHeaderPtr(
        sam_hdr_read(htsFilePtr.get()), bam_hdr_destroy);
    if (!htsHeaderPtr)
    {
        throw std::runtime_error("Failed to read the header of " + htsFilePath);
    }

    return htshelpers::decodeContigInfo(htsHeaderPtr.get());
}

}
