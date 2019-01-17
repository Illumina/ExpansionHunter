//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>
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

#include <memory>
#include <string>
#include <unordered_map>

// cppcheck-suppress missingInclude
#include "htslib/hts.h"
// cppcheck-suppress missingInclude
#include "htslib/sam.h"

#include "graphalign/GraphAlignment.hh"
#include "graphcore/GraphReferenceMapping.hh"
#include "graphio/AlignmentWriter.hh"

#include "common/ReferenceContigInfo.hh"
#include "reads/Read.hh"
#include "region_spec/LocusSpecification.hh"

namespace ehunter
{

class BamletWriter : public graphtools::AlignmentWriter
{
public:
    BamletWriter(
        const std::string& bamletPath, const ReferenceContigInfo& contigInfo, const RegionCatalog& regionCatalog);
    ~BamletWriter() override = default;

    void write(
        const std::string& locusId, const std::string& fragmentName, const std::string& query, bool isFirstMate,
        const graphtools::GraphAlignment& alignment) override;

private:
    void writeHeader();
    void write(
        const graphtools::ReferenceInterval& interval, const std::string& fragmentName, const std::string& query,
        bool isFirstMate, const graphtools::GraphAlignment& alignment);

    std::unique_ptr<htsFile, decltype(&hts_close)> filePtr_;
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> bamHeader_;
    ReferenceContigInfo contigInfo_;

    std::unordered_map<std::string, graphtools::GraphReferenceMapping> graphReferenceMappings_;
};

}
