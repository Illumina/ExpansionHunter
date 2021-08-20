//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Include the fai class from samtools
#include "htslib/faidx.h"

#include "core/GenomicRegion.hh"
#include "core/ReferenceContigInfo.hh"

namespace ehunter
{

class Reference
{
public:
    /**
     * @param contigName Name of the reference contig
     * @param start 0-based, inclusive
     * @param end 0-based, exclusive
     * @return Reference sequence in upper case
     */
    virtual std::string getSequence(const std::string& contigName, int64_t start, int64_t end) const = 0;

    virtual std::string getSequence(const GenomicRegion& region) const = 0;

    virtual const ReferenceContigInfo& contigInfo() const = 0;
};

/**
 * Reference genome implementation backed by a FASTA file read through HTSlib
 */
class FastaReference : public Reference
{
public:
    explicit FastaReference(const std::string& referencePath, const ReferenceContigInfo& contigInfo);
    ~FastaReference();

    std::string getSequence(const std::string& contigIndex, int64_t start, int64_t end) const override;
    std::string getSequence(const GenomicRegion& region) const override;

    const ReferenceContigInfo& contigInfo() const override { return bamHeaderContigInfo_; }

private:
    // A stub for verifying consistency of contig information in FASTA index and BAM header
    void assertConsistency() const;

    std::string referencePath_;
    faidx_t* htsFastaIndexPtr_;

    // fastaContigInfo_ is meant for internal use only; it should not be exposed though public interface of the class
    ReferenceContigInfo fastaContigInfo_;
    // bamHeaderContigInfo_ is meant to be exposed through the public interface
    ReferenceContigInfo bamHeaderContigInfo_;
};

}
