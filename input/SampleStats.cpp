//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "input/SampleStats.hh"

#include <boost/filesystem.hpp>

#include "htslib/hts.h"
#include "htslib/sam.h"

using std::string;

namespace ehunter
{

int extractReadLength(const string& bamPath)
{
    // Open a BAM file for reading.
    samFile* htsFilePtr = sam_open(bamPath.c_str(), "r");
    if (!htsFilePtr)
    {
        throw std::runtime_error("Failed to read BAM file '" + bamPath + "'");
    }
    bam_hdr_t* htsHeaderPtr = sam_hdr_read(htsFilePtr);
    if (!htsHeaderPtr)
    {
        throw std::runtime_error("BamFile::Init: Failed to read BAM header: '" + bamPath + "'");
    }

    enum
    {
        kSupplementaryAlign = 0x800,
        kSecondaryAlign = 0x100
    };

    int readLength = 99;
    bam1_t* htsAlignmentPtr = bam_init1();
    int ret;
    while ((ret = sam_read1(htsFilePtr, htsHeaderPtr, htsAlignmentPtr)) >= 0)
    {
        const bool isSupplementary = htsAlignmentPtr->core.flag & kSupplementaryAlign;
        const bool isSecondary = htsAlignmentPtr->core.flag & kSecondaryAlign;
        const bool isPrimaryAlign = (!isSupplementary) && (!isSecondary);
        if (isPrimaryAlign)
        {
            readLength = htsAlignmentPtr->core.l_qseq;
            break;
        }
    }

    if (ret < 0)
    {
        throw std::runtime_error("Failed to extract a read from BAM file");
    }

    bam_destroy1(htsAlignmentPtr);
    bam_hdr_destroy(htsHeaderPtr);
    sam_close(htsFilePtr);

    return readLength;
}

bool isBamFile(const string& htsFilePath)
{
    string extension = boost::filesystem::extension(htsFilePath);
    return extension == ".bam";
}

}
