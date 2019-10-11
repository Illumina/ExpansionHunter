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

#include "common/GenomicRegion.hh"
#include "common/ReferenceContigInfo.hh"
#include "reads/Read.hh"

namespace ehunter
{

namespace htshelpers
{

    class HtsFileSeeker
    {
    public:
        HtsFileSeeker(const std::string& htsFilePath, const std::string& htsReferencePath);
        ~HtsFileSeeker();
        void setRegion(const GenomicRegion& region);
        bool trySeekingToNextPrimaryAlignment();

        int currentReadChromIndex() const;
        const std::string& currentReadChrom() const;
        int64_t currentReadPosition() const;
        int currentMateChromIndex() const;
        const std::string& currentMateChrom() const;
        int64_t currentMatePosition() const;

        MappedRead decodeRead() const;

    private:
        enum class Status
        {
            kStreamingReads,
            kFinishedStreaming
        };

        void openFile();
        void loadHeader();
        void loadIndex();
        void closeRegion();

        const std::string htsFilePath_;
        const std::string htsReferencePath_;
        ReferenceContigInfo contigInfo_;
        Status status_ = Status::kFinishedStreaming;

        htsFile* htsFilePtr_ = nullptr;
        bam_hdr_t* htsHeaderPtr_ = nullptr;
        hts_idx_t* htsIndexPtr_ = nullptr;
        hts_itr_t* htsRegionPtr_ = nullptr;
        bam1_t* htsAlignmentPtr_ = nullptr;
    };

}

}
