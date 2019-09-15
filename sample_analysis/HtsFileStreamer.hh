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

#include "common/ReferenceContigInfo.hh"
#include "reads/Read.hh"

namespace ehunter
{

namespace htshelpers
{

    class HtsFileStreamer
    {
    public:
        HtsFileStreamer(const std::string& htsFilePath)
            : htsFilePath_(htsFilePath)
            , contigInfo_({})
        {
            openHtsFile();
            loadHeader();
            prepareForStreamingAlignments();
        }
        ~HtsFileStreamer();

        bool trySeekingToNextPrimaryAlignment();

        int32_t currentReadContigId() const;
        int32_t currentReadPosition() const;
        int32_t currentReadLength() const;
        int32_t currentMateContigId() const;
        int32_t currentMatePosition() const;

        bool isStreamingAlignedReads() const;

        MappedRead decodeRead() const;

    private:
        enum class Status
        {
            kStreamingReads,
            kFinishedStreaming
        };

        void openHtsFile();
        void loadHeader();
        void prepareForStreamingAlignments();

        const std::string htsFilePath_;
        ReferenceContigInfo contigInfo_;
        Status status_ = Status::kStreamingReads;

        htsFile* htsFilePtr_ = nullptr;
        bam1_t* htsAlignmentPtr_ = nullptr;
        bam_hdr_t* htsHeaderPtr_ = nullptr;
    };

}
}
