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

#include "sample_analysis/HtsFileSeeker.hh"

#include <cstdint>
#include <stdexcept>

#include "sample_analysis/HtsHelpers.hh"

using std::string;

namespace ehunter
{

namespace htshelpers
{

    HtsFileSeeker::HtsFileSeeker(const string& htsFilePath)
        : htsFilePath_(htsFilePath)
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
    }

    void HtsFileSeeker::loadHeader()
    {
        htsHeaderPtr_ = sam_hdr_read(htsFilePtr_);

        if (!htsHeaderPtr_)
        {
            throw std::runtime_error("Failed to read header of " + htsFilePath_);
        }

        const int32_t numContigs = htsHeaderPtr_->n_targets;

        for (int32_t contigInd = 0; contigInd != numContigs; ++contigInd)
        {
            const string contig = htsHeaderPtr_->target_name[contigInd];
            contigNames_.push_back(contig);
        }
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

    void HtsFileSeeker::setRegion(const Region& region)
    {
        closeRegion();

        const string regionEncoding = region.ToString();
        htsRegionPtr_ = sam_itr_querys(htsIndexPtr_, htsHeaderPtr_, regionEncoding.c_str());

        if (htsRegionPtr_ == nullptr)
        {
            throw std::runtime_error("Failed to extract reads from " + regionEncoding);
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

    reads::Read HtsFileSeeker::decodeRead(reads::LinearAlignmentStats& alignmentStats) const
    {
        reads::Read read;
        DecodeAlignedRead(htsAlignmentPtr_, read, alignmentStats);

        return read;
    }

}

}
