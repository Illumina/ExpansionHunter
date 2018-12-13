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
    class MateExtractor
    {
    public:
        MateExtractor(const std::string& htsFilePath);
        ~MateExtractor();

        reads::Read extractMate(const reads::Read& read, const reads::LinearAlignmentStats& alignmentStats);

    private:
        void openFile();
        void loadHeader();
        void loadIndex();

        const std::string htsFilePath_;
        std::vector<std::string> contigNames_;

        htsFile* htsFilePtr_ = nullptr;
        bam_hdr_t* htsHeaderPtr_ = nullptr;
        hts_idx_t* htsIndexPtr_ = nullptr;
        bam1_t* htsAlignmentPtr_ = nullptr;
    };

}

}
