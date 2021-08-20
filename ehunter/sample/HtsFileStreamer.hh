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

#include <string>
#include <vector>

extern "C"
{
#include "htslib/hts.h"
#include "htslib/sam.h"
}

#include "core/Read.hh"
#include "core/ReferenceContigInfo.hh"

namespace ehunter
{

namespace htshelpers
{

class HtsFileStreamer
{
public:
    /// \param[in] decompressionThreads Total size of thread pool used for bgzip decompression of the hts file. Has no
    /// effect if the file is uncompressed. When set to one or less the calling thread handles all decompression and no
    /// thread pool is used.
    ///
    HtsFileStreamer(
        const std::string& htsFilePath, const std::string& htsReferencePath, const unsigned decompressionThreads = 1)
        : htsFilePath_(htsFilePath)
        , htsReferencePath_(htsReferencePath)
        , contigInfo_({})
    {
        openHtsFile(decompressionThreads);
        loadHeader();
        prepareForStreamingAlignments();
    }
    ~HtsFileStreamer();

    bool trySeekingToNextPrimaryAlignment();

    int32_t currentReadContigId() const { return htsAlignmentPtr_->core.tid; }
    hts_pos_t currentReadPosition() const { return htsAlignmentPtr_->core.pos; }
    int32_t currentReadLength() const;
    int32_t currentMateContigId() const { return htsAlignmentPtr_->core.mtid; }
    hts_pos_t currentMatePosition() const { return htsAlignmentPtr_->core.mpos; }

    bool currentIsPaired() const { return (htsAlignmentPtr_->core.flag & BAM_FPAIRED); }

    bool isStreamingAlignedReads() const;

    Read decodeRead() const;

private:
    enum class Status
    {
        kStreamingReads,
        kFinishedStreaming
    };

    void openHtsFile(unsigned decompressionThreads);
    void loadHeader();
    void prepareForStreamingAlignments();

    const std::string htsFilePath_;
    const std::string htsReferencePath_;
    ReferenceContigInfo contigInfo_;
    Status status_ = Status::kStreamingReads;

    htsFile* htsFilePtr_ = nullptr;
    bam1_t* htsAlignmentPtr_ = nullptr;
    bam_hdr_t* htsHeaderPtr_ = nullptr;
    htsThreadPool htsThreadPool_ = { nullptr, 0 };
};

}
}
