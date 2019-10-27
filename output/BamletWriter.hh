//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>
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
#include "locus_spec/GraphLocusSpecification.hh"
#include "locus_spec/LocusSpecification.hh"
#include "reads/Read.hh"

namespace ehunter
{

class BamletWriter : public graphtools::AlignmentWriter
{
public:
    BamletWriter(
        const std::string& bamletPath, const ReferenceContigInfo& contigInfo, const LocusCatalog& regionCatalog);
    ~BamletWriter() override = default;

    void write(
        const std::string& locusId, const std::string& fragmentName, const std::string& query, bool isFirstMate,
        bool isReversed, bool isMateReversed, const graphtools::GraphAlignment& alignment) override;

private:
    void writeHeader();
    void write(
        const graphtools::ReferenceInterval& interval, const std::string& fragmentName, const std::string& query,
        bool isFirstMate, bool isReversed, bool isMateReversed, const graphtools::GraphAlignment& alignment);

    std::unique_ptr<htsFile, decltype(&hts_close)> filePtr_;
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> bamHeader_;
    ReferenceContigInfo contigInfo_;

    std::unordered_map<std::string, graphtools::GraphReferenceMapping> graphReferenceMappings_;
};

using BamletWriterPtr = std::shared_ptr<BamletWriter>;
}
