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

// Include the fai class from samtools
#include "htslib/faidx.h"

#include "common/GenomicRegion.hh"

namespace ehunter
{

using pos_t = size_t;

class Reference
{
public:
    /**
     * @param chrom Name of the reference contig (chromosome)
     * @param start 0-based, inclusive
     * @param end 0-based, exclusive
     * @return Reference sequence in upper case
     */
    virtual std::string getSequence(const std::string& chrom, pos_t start, pos_t end) const = 0;
};

/**
 * Reference Genome implementation backed by a fasta file read through htslib
 */
class FastaReference : public Reference
{
public:
    explicit FastaReference(const std::string& genome_path);
    ~FastaReference();

    std::string getSequence(const std::string& chrom, pos_t start, pos_t end) const override;

private:
    std::string genome_path_;
    faidx_t* fai_ptr_;
};

}
