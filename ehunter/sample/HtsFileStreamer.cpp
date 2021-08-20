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

#include "sample/HtsFileStreamer.hh"

#include "core/HtsHelpers.hh"

using std::string;

namespace ehunter
{

namespace htshelpers
{

void HtsFileStreamer::openHtsFile(const unsigned decompressionThreads)
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

    /// Create thread pool for bgzf block decompression. The behavior of htslib seems to be to use this pool instead of
    /// (not in addition to) the calling thread, therefore there is no point in creating a decompression thread-pool
    /// with less than 2 threads.
    if (decompressionThreads > 1)
    {
        htsThreadPool_.pool = hts_tpool_init(decompressionThreads);
        if (not htsThreadPool_.pool)
        {
            throw std::runtime_error(
                "HtsFileStreamer: Error creating htslib threadpool with " + std::to_string(decompressionThreads)
                + " threads.");
        }
        hts_set_opt(htsFilePtr_, HTS_OPT_THREAD_POOL, &htsThreadPool_);
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

bool HtsFileStreamer::isStreamingAlignedReads() const
{
    return status_ != Status::kFinishedStreaming && currentReadContigId() != -1;
}

Read HtsFileStreamer::decodeRead() const { return htshelpers::decodeRead(htsAlignmentPtr_); }

HtsFileStreamer::~HtsFileStreamer()
{
    bam_destroy1(htsAlignmentPtr_);
    htsAlignmentPtr_ = nullptr;

    bam_hdr_destroy(htsHeaderPtr_);
    htsHeaderPtr_ = nullptr;

    sam_close(htsFilePtr_);
    htsFilePtr_ = nullptr;

    if (htsThreadPool_.pool)
    {
        hts_tpool_destroy(htsThreadPool_.pool);
    }
}

}

}
