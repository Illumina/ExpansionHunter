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

#include <boost/optional.hpp>

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
class MateExtractor
{
public:
    MateExtractor(const std::string& htsFilePath, const std::string& htsReferencePath);
    ~MateExtractor();

    boost::optional<Read>
    extractMate(const Read& read, const LinearAlignmentStats& alignmentStats, LinearAlignmentStats& mateStats);

private:
    void openFile();
    void loadHeader();
    void loadIndex();

    std::string htsFilePath_;
    std::string htsReferencePath_;
    ReferenceContigInfo contigInfo_;

    htsFile* htsFilePtr_ = nullptr;
    bam_hdr_t* htsHeaderPtr_ = nullptr;
    hts_idx_t* htsIndexPtr_ = nullptr;
    bam1_t* htsAlignmentPtr_ = nullptr;
};

}

}
