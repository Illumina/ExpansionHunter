//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
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

#pragma once

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// Include the fai class from samtools
#include "htslib/faidx.h"

#include "common/GenomicRegion.hh"
#include "common/ReferenceContigInfo.hh"

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

    const ReferenceContigInfo& contigInfo() const override { return contigInfo_; }

private:
    std::string referencePath_;
    faidx_t* htsFastaIndexPtr_;

    ReferenceContigInfo contigInfo_;
};

}
