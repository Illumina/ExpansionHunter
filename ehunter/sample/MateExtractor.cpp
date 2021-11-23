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

#include "sample/MateExtractor.hh"

#include <stdexcept>

#include "core/HtsHelpers.hh"

namespace ehunter
{

using boost::optional;
using std::string;

namespace htshelpers
{
MateExtractor::MateExtractor(const string& htsFilePath, const std::string& htsReferencePath)
    : htsFilePath_(htsFilePath)
    , htsReferencePath_(htsReferencePath)
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

    // Required step for parsing of some CRAMs
    if (hts_set_fai_filename(htsFilePtr_, htsReferencePath_.c_str()) != 0)
    {
        throw std::runtime_error("Failed to set index of: " + htsReferencePath_);
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

std::vector<std::pair<Read, LinearAlignmentStats>> MateExtractor::extractMates(
    const MateRegionToRecover& mateRegionToRecover)
{

    std::unordered_set<ReadId, boost::hash<ReadId>> mateReadIdsNotFoundYet = mateRegionToRecover.mateReadIds;
    std::vector<std::pair<Read, LinearAlignmentStats>> extractedMatesAndAlignmentStats;

    const int32_t searchRegionContigIndex = mateRegionToRecover.genomicRegion.contigIndex();
    const int32_t searchRegionStart = mateRegionToRecover.genomicRegion.start();
    const int32_t searchRegionEnd = mateRegionToRecover.genomicRegion.end();

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
        const bool isSecondaryAlignment = htsAlignmentPtr_->core.flag & BAM_FSECONDARY;
        const bool isSupplementaryAlignment = htsAlignmentPtr_->core.flag & BAM_FSUPPLEMENTARY;
        const bool isPrimaryAlignment = !(isSecondaryAlignment || isSupplementaryAlignment);
        if (!isPrimaryAlignment)
        {
            continue;
        }

        ReadId putativeMateReadId = decodeReadId(htsAlignmentPtr_);

        // if this is one of the fragment ids we're looking for, check if forms a proper pair with the 1st read with
        // this fragmentId. If it does, decode this read and add it as a mate
        bool foundRequestedMate = mateReadIdsNotFoundYet.count(putativeMateReadId) > 0;
        if (foundRequestedMate) {
            LinearAlignmentStats mateStats = decodeAlignmentStats(htsAlignmentPtr_);
            Read putativeMate = htshelpers::decodeRead(htsAlignmentPtr_);
            auto mateAndAlignmentStats = std::make_pair(putativeMate, mateStats);
            extractedMatesAndAlignmentStats.push_back(mateAndAlignmentStats);
            mateReadIdsNotFoundYet.erase(putativeMate.readId());
            if (mateReadIdsNotFoundYet.size() == 0) {
                break;  // found all requested mates
            }
        }
    }
    hts_itr_destroy(htsRegionPtr_);

    if (mateReadIdsNotFoundYet.size() > 0) {
        throw std::runtime_error("Failed to recover " + std::to_string(mateReadIdsNotFoundYet.size()) + " mates");
    }

    return extractedMatesAndAlignmentStats;
}

}

}
