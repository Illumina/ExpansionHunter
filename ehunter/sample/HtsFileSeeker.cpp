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

#include "sample/HtsFileSeeker.hh"

#include <cstdint>
#include <stdexcept>

#include "core/HtsHelpers.hh"

using std::string;

namespace ehunter
{

namespace htshelpers
{

HtsFileSeeker::HtsFileSeeker(const string& htsFilePath, const std::string& htsReferencePath)
    : htsFilePath_(htsFilePath)
    , htsReferencePath_(htsReferencePath)
    , contigInfo_({})
{
    openFile();
    loadHeader();
    loadIndex();
    htsAlignmentPtr_ = bam_init1();
}

HtsFileSeeker::~HtsFileSeeker()
{
    bam_destroy1(htsAlignmentPtr_);
    htsAlignmentPtr_ = nullptr;

    if (htsRegionPtr_)
    {
        hts_itr_destroy(htsRegionPtr_);
        htsRegionPtr_ = nullptr;
    }

    hts_idx_destroy(htsIndexPtr_);
    htsIndexPtr_ = nullptr;

    bam_hdr_destroy(htsHeaderPtr_);
    htsHeaderPtr_ = nullptr;

    sam_close(htsFilePtr_);
    htsFilePtr_ = nullptr;
}

void HtsFileSeeker::openFile()
{
    htsFilePtr_ = sam_open(htsFilePath_.c_str(), "r");

    if (!htsFilePtr_)
    {
        throw std::runtime_error("Failed to read BAM file " + htsFilePath_);
    }

    // Required step for parsing of some CRAMs
    if (hts_set_fai_filename(htsFilePtr_, htsReferencePath_.c_str()) != 0)
    {
        throw std::runtime_error("Failed to set index of: " + htsReferencePath_);
    }
}

void HtsFileSeeker::loadHeader()
{
    htsHeaderPtr_ = sam_hdr_read(htsFilePtr_);

    if (!htsHeaderPtr_)
    {
        throw std::runtime_error("Failed to read header of " + htsFilePath_);
    }

    contigInfo_ = htshelpers::decodeContigInfo(htsHeaderPtr_);
}

void HtsFileSeeker::loadIndex()
{
    htsIndexPtr_ = sam_index_load(htsFilePtr_, htsFilePath_.c_str());

    if (!htsIndexPtr_)
    {
        throw std::runtime_error("Failed to read index of " + htsFilePath_);
    }
}

void HtsFileSeeker::closeRegion()
{
    if (htsRegionPtr_)
    {
        hts_itr_destroy(htsRegionPtr_);
        htsRegionPtr_ = nullptr;
    }
}

void HtsFileSeeker::setRegion(const GenomicRegion& region)
{
    closeRegion();

    htsRegionPtr_ = sam_itr_queryi(htsIndexPtr_, region.contigIndex(), region.start(), region.end());

    if (htsRegionPtr_ == nullptr)
    {
        throw std::runtime_error("Failed to extract reads from " + encode(contigInfo_, region));
    }

    status_ = Status::kStreamingReads;
}

bool HtsFileSeeker::trySeekingToNextPrimaryAlignment()
{
    if (status_ != Status::kStreamingReads)
    {
        return false;
    }

    int32_t returnCode = 0;

    while ((returnCode = sam_itr_next(htsFilePtr_, htsRegionPtr_, htsAlignmentPtr_)) >= 0)
    {
        if (isPrimaryAlignment(htsAlignmentPtr_))
            return true;
    }

    status_ = Status::kFinishedStreaming;

    if (returnCode < -1)
    {
        throw std::runtime_error("Failed to extract a record from " + htsFilePath_);
    }

    return false;
}

Read HtsFileSeeker::decodeRead(LinearAlignmentStats& alignmentStats) const
{
    alignmentStats = decodeAlignmentStats(htsAlignmentPtr_);
    return htshelpers::decodeRead(htsAlignmentPtr_);
}

}

}
