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

#include "sample_analysis/MateExtractor.hh"

#include <stdexcept>

#include "common/HtsHelpers.hh"

namespace ehunter
{

using boost::optional;
using std::string;

namespace htshelpers
{
    MateExtractor::MateExtractor(const string& htsFilePath)
        : htsFilePath_(htsFilePath)
        , contigInfo_({})
    {
        openFile();
        loadHeader();
        loadIndex();
        htsAlignmentPtr_ = bam_init1();
    }

    MateExtractor::~MateExtractor()
    {
        bam_destroy1(htsAlignmentPtr_);
        htsAlignmentPtr_ = nullptr;

        hts_idx_destroy(htsIndexPtr_);
        htsIndexPtr_ = nullptr;

        bam_hdr_destroy(htsHeaderPtr_);
        htsHeaderPtr_ = nullptr;

        sam_close(htsFilePtr_);
        htsFilePtr_ = nullptr;
    }

    void MateExtractor::openFile()
    {
        htsFilePtr_ = sam_open(htsFilePath_.c_str(), "r");

        if (!htsFilePtr_)
        {
            throw std::runtime_error("Failed to read BAM file " + htsFilePath_);
        }
    }

    void MateExtractor::loadHeader()
    {
        htsHeaderPtr_ = sam_hdr_read(htsFilePtr_);

        if (!htsHeaderPtr_)
        {
            throw std::runtime_error("Failed to read header of " + htsFilePath_);
        }

        contigInfo_ = htshelpers::decodeContigInfo(htsHeaderPtr_);
    }

    void MateExtractor::loadIndex()
    {
        htsIndexPtr_ = sam_index_load(htsFilePtr_, htsFilePath_.c_str());

        if (!htsIndexPtr_)
        {
            throw std::runtime_error("Failed to read index of " + htsFilePath_);
        }
    }

    optional<Read> MateExtractor::extractMate(const Read& read, const LinearAlignmentStats& alignmentStats)
    {
        const int32_t searchRegionContigIndex
            = alignmentStats.isMateMapped ? alignmentStats.mateChromId : alignmentStats.chromId;

        const int32_t searchRegionStart = alignmentStats.isMateMapped ? alignmentStats.matePos : alignmentStats.pos;
        const int32_t searchRegionEnd = searchRegionStart + 1;

        hts_itr_t* htsRegionPtr_
            = sam_itr_queryi(htsIndexPtr_, searchRegionContigIndex, searchRegionStart, searchRegionEnd);

        if (!htsRegionPtr_)
        {
            const string& contigName = contigInfo_.getContigName(searchRegionContigIndex);
            const string regionEncoding
                = contigName + ":" + std::to_string(searchRegionStart) + "-" + std::to_string(searchRegionEnd);

            throw std::logic_error("Unable to jump to " + regionEncoding + " to recover a mate");
        }

        while (sam_itr_next(htsFilePtr_, htsRegionPtr_, htsAlignmentPtr_) >= 0)
        {
            Read putativeMate = htshelpers::decodeRead(htsAlignmentPtr_);

            const bool belongToSameFragment = read.fragmentId() == putativeMate.fragmentId();
            const bool formProperPair = read.mateNumber() != putativeMate.mateNumber();
            if (belongToSameFragment && formProperPair)
            {
                hts_itr_destroy(htsRegionPtr_);
                return putativeMate;
            }
        }
        hts_itr_destroy(htsRegionPtr_);

        return optional<Read>();
    }

}

}
