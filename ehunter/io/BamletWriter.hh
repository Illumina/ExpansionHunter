//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>,
//         Chris Saunders <csaunders@illumina.com>
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
#include <thread>
#include <unordered_map>

#include "boost/noncopyable.hpp"
// cppcheck-suppress missingInclude
#include "htslib/hts.h"
// cppcheck-suppress missingInclude
#include "htslib/sam.h"

#include "core/ConcurrentQueue.hh"
#include "core/Read.hh"
#include "core/ReferenceContigInfo.hh"
#include "graphalign/GraphAlignment.hh"
#include "graphcore/GraphReferenceMapping.hh"
#include "graphio/AlignmentWriter.hh"
#include "locus/LocusSpecification.hh"

namespace ehunter
{

/// Supports multiple threads calling the write() method on the same object. To make this more efficient for high
/// thread counts, this object creates its own asynchronous write thread to prevent calling threads from blocking
/// on the final bam record write operation.
///
class BamletWriter : public graphtools::AlignmentWriter, private boost::noncopyable
{
public:
    BamletWriter(
        const std::string& bamletPath, const ReferenceContigInfo& contigInfo, const RegionCatalog& regionCatalog);
    ~BamletWriter() override;

    /// Thread safe
    void write(
        const std::string& locusId, const std::string& fragmentName, const std::string& query, bool isFirstMate,
        bool isReversed, bool isMateReversed, const graphtools::GraphAlignment& alignment) override;

private:
    void writeHeader();

    /// Thread safe
    void write(
        const graphtools::ReferenceInterval& interval, const std::string& fragmentName, const std::string& query,
        bool isFirstMate, bool isReversed, bool isMateReversed, const graphtools::GraphAlignment& alignment);

    /// Function executed by dedicated bam writer thread
    void writeHtsAlignments();

    std::unique_ptr<htsFile, decltype(&hts_close)> filePtr_;
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> bamHeader_;
    ReferenceContigInfo contigInfo_;

    std::unordered_map<std::string, graphtools::GraphReferenceMapping> graphReferenceMappings_;

    ConcurrentQueue<bam1_t*> writeQueue_;
    std::thread writeThread_;
};

}
