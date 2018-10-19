// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include <memory>
#include <string>
#include <vector>

// cppcheck-suppress missingInclude
#include "htslib/hts.h"
// cppcheck-suppress missingInclude
#include "htslib/sam.h"

#include "graphalign/GraphAlignment.hh"
#include "graphcore/GraphReferenceMapping.hh"

using namespace graphtools;
using Sequence = std::string; // Type of read sequences
namespace graphIO
{

using ReferenceContigs = std::vector<std::pair<std::string, uint32_t>>;

/**
 * Subset of information on a BAM record for graph-alignment output
 */
struct BamAlignment
{
    std::string chromName; // Has to match a contig in the BAM header
    int pos = -1; // 0-based
    bool isPaired = false;
    bool isMate1 = false;
    std::string fragmentName;
    Sequence sequence;
    std::vector<int> BaseQualities;
    std::string graphCigar; // Represents the graph alignment of the read (in a string BAM tag)
};

/**
 * Write Graph-alignments to a BAM file
 */
class BamWriter
{
public:
    /**
     * Paired-end status of an alignment
     */
    enum class PairingInfo
    {
        Unpaired,
        FirstMate,
        SecondMate
    };

    /**
     * @param bamPath Path to output BAM file
     * @param contigs Contig names and sequence lengths for BAM header
     * @throws If BamFile cannot be created
     */
    explicit BamWriter(std::string const& bamPath, ReferenceContigs& contigs);
    /**
     * Create an unplaced BAM alignment with a graph CIGAR tag
     * @param fragmentName Name of the aligned read (fragment)
     * @param sequence Sequence of the aligned read
     * @param qualities BaseQ. Must be empty or have the same length as sequence
     * @param pairing If paired-end, which mate is this?
     * @param graphAlign Graph-CIGAR string of alignment
     */
    BamAlignment makeAlignment(
        std::string const& fragmentName, Sequence const& sequence, std::vector<int> const& qualities,
        BamWriter::PairingInfo pairing, std::string const& graphAlign) const;

    /**
     * Project a graph alignment to the reference genome and output as placed but unmapped BAM record
     * @param refMap Graph to refernece genome coordinate projection
     * @param fragmentName Name of the aligned read (fragment)
     * @param sequence Sequence of the aligned read
     * @param qualities BaseQ. Must be empty or have the same length as sequence
     * @param pairing If paired-end, which mate is this?
     * @param graphAlign Graph-aligment to project and output
     */
    BamAlignment makeAlignment(
        GraphReferenceMapping const& refMap, std::string const& fragmentName, Sequence const& sequence,
        std::vector<int> const& qualities, BamWriter::PairingInfo pairing, GraphAlignment const& align) const;
    /**
     * Write a BAM alignment
     * @param align Alignment to write
     */
    void writeAlignment(BamAlignment& align);

private:
    int writeHeader(std::string const& initHeader, ReferenceContigs& contigs);

    std::unique_ptr<htsFile, decltype(&hts_close)> filePtr_;
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> bamHeader_;

    static std::string const initHeader; // Dummy header line
    static std::string const graphCigarBamTag; // Custom tag to use for graphCIGAR string
};
}