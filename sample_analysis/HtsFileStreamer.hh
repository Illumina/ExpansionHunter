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

#pragma once

#include <string>
#include <vector>

extern "C"
{
#include "htslib/hts.h"
#include "htslib/sam.h"
}

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
        {
            openHtsFile();
            loadHeader();
            prepareForStreamingAlignments();
        }
        ~HtsFileStreamer();

        bool trySeekingToNextPrimaryAlignment();

        int32_t currentReadChromIndex() const;
        const std::string& currentReadChrom() const;
        int32_t currentReadPosition() const;
        int32_t currentMateChromIndex() const;
        const std::string& currentMateChrom() const;
        int32_t currentMatePosition() const;

        bool isStreamingAlignedReads() const;

        reads::Read decodeRead() const;

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
        std::vector<std::string> chromNames_;
        Status status_ = Status::kStreamingReads;

        htsFile* htsFilePtr_ = nullptr;
        bam1_t* htsAlignmentPtr_ = nullptr;
        bam_hdr_t* htsHeaderPtr_ = nullptr;
    };

}
}
