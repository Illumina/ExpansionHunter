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

#include "sample_analysis/HtsFileStreamer.hh"

#include "common/HtsHelpers.hh"

using std::string;

namespace ehunter
{

namespace htshelpers
{

    void HtsFileStreamer::openHtsFile()
    {
        htsFilePtr_ = sam_open(htsFilePath_.c_str(), "r");

        if (!htsFilePtr_)
        {
            throw std::runtime_error("Failed to read BAM file " + htsFilePath_);
        }
    }

    void HtsFileStreamer::loadHeader()
    {
        htsHeaderPtr_ = sam_hdr_read(htsFilePtr_);

        if (!htsHeaderPtr_)
        {
            throw std::runtime_error("Failed to read header of " + htsFilePath_);
        }

        contigInfo_ = htshelpers::decodeContigInfo(htsHeaderPtr_);
    }

    void HtsFileStreamer::prepareForStreamingAlignments() { htsAlignmentPtr_ = bam_init1(); }

    bool HtsFileStreamer::trySeekingToNextPrimaryAlignment()
    {
        if (status_ != Status::kStreamingReads)
        {
            return false;
        }

        int32_t returnCode = 0;

        while ((returnCode = sam_read1(htsFilePtr_, htsHeaderPtr_, htsAlignmentPtr_)) >= 0)
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

    int32_t HtsFileStreamer::currentReadContigId() const { return htsAlignmentPtr_->core.tid; }
    int32_t HtsFileStreamer::currentReadPosition() const { return htsAlignmentPtr_->core.pos; }
    int32_t HtsFileStreamer::currentMateContigId() const { return htsAlignmentPtr_->core.mtid; }
    int32_t HtsFileStreamer::currentMatePosition() const { return htsAlignmentPtr_->core.mpos; }

    bool HtsFileStreamer::isStreamingAlignedReads() const
    {
        return status_ != Status::kFinishedStreaming && currentReadContigId() != -1;
    }

    MappedRead HtsFileStreamer::decodeRead() const { return htshelpers::decodeRead(htsAlignmentPtr_); }

    HtsFileStreamer::~HtsFileStreamer()
    {
        bam_destroy1(htsAlignmentPtr_);
        htsAlignmentPtr_ = nullptr;

        bam_hdr_destroy(htsHeaderPtr_);
        htsHeaderPtr_ = nullptr;

        sam_close(htsFilePtr_);
        htsFilePtr_ = nullptr;
    }

}

}
