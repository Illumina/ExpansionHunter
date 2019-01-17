//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
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

#include "sample_analysis/HtsFileStreamer.hh"

#include "sample_analysis/HtsHelpers.hh"

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

        const int32_t numChroms = htsHeaderPtr_->n_targets;

        for (int32_t chromInd = 0; chromInd != numChroms; ++chromInd)
        {
            const string chrom = htsHeaderPtr_->target_name[chromInd];
            chromNames_.push_back(chrom);
        }
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
            const bool isPrimaryAlignment = !(htsAlignmentPtr_->core.flag & htshelpers::kIsNotPrimaryLine);

            if (isPrimaryAlignment)
            {
                return true;
            }
        }

        status_ = Status::kFinishedStreaming;

        if (returnCode < -1)
        {
            throw std::runtime_error("Failed to extract a record from " + htsFilePath_);
        }

        return false;
    }

    int32_t HtsFileStreamer::currentReadChromIndex() const { return htsAlignmentPtr_->core.tid; }
    const std::string& HtsFileStreamer::currentReadChrom() const { return chromNames_[currentReadChromIndex()]; }
    int32_t HtsFileStreamer::currentReadPosition() const { return htsAlignmentPtr_->core.pos; }

    int32_t HtsFileStreamer::currentMateChromIndex() const { return htsAlignmentPtr_->core.mtid; }
    const std::string& HtsFileStreamer::currentMateChrom() const { return chromNames_[currentMateChromIndex()]; }
    int32_t HtsFileStreamer::currentMatePosition() const { return htsAlignmentPtr_->core.mpos; }

    bool HtsFileStreamer::isStreamingAlignedReads() const
    {
        return status_ != Status::kFinishedStreaming && currentReadChromIndex() != -1;
    }

    reads::Read HtsFileStreamer::decodeRead() const
    {
        reads::Read read;
        reads::LinearAlignmentStats alignmentStats;
        DecodeAlignedRead(htsAlignmentPtr_, read, alignmentStats);

        return read;
    }

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
