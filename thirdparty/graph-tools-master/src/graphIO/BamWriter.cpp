//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>
//
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

#include "graphIO/BamWriter.hh"

#include <iostream>
#include <sstream>
#include <stdexcept>

// cppcheck-suppress missingInclude
#include "htslib/bgzf.h"
// cppcheck-suppress missingInclude
#include "htslib/kseq.h"
// cppcheck-suppress missingInclude
#include "htslib/sam.h"

using std::string;

// Set a base in the bit-backed read sequence representation of a BAM record. (4 bits per base).
// Sets position i to base c in sequence s.
#define bam1_seq_seti(s, i, c) ((s)[(i) >> 1] = ((s)[(i) >> 1] & 0xf << (((i)&1) << 2)) | (c) << ((~(i)&1) << 2))

namespace graphIO
{

string const BamWriter::initHeader = "@HD\tVN:1.4\tSO:unknown\n";
string const BamWriter::graphCigarBamTag = "XG";

BamWriter::BamWriter(string const& bamPath, ReferenceContigs& contigs)
    : filePtr_(hts_open(bamPath.c_str(), "wb"), hts_close)
    , bamHeader_(bam_hdr_init(), bam_hdr_destroy)
{
    if (writeHeader(initHeader, contigs) != 0)
    {
        throw std::logic_error("Failed to write header");
    }
}

int BamWriter::writeHeader(std::string const& initHeader, ReferenceContigs& contigs)
{
    bamHeader_->l_text = strlen(initHeader.c_str());
    bamHeader_->text = strdup(initHeader.c_str());
    bamHeader_->n_targets = contigs.size();

    // All this memory gets freed by the header (bam_hdr_destroy)
    bamHeader_->target_name = new char*[contigs.size()];
    bamHeader_->target_len = new uint32_t[contigs.size()];
    for (size_t i = 0; i != contigs.size(); ++i)
    {
        bamHeader_->target_name[i] = new char[contigs[i].first.length() + 1];
        memcpy(bamHeader_->target_name[i], contigs[i].first.c_str(), contigs[i].first.length() + 1);
        bamHeader_->target_len[i] = contigs[i].second;
    }
    return bam_hdr_write(filePtr_->fp.bgzf, bamHeader_.get());
}

void BamWriter::writeAlignment(BamAlignment& align)
{
    bam1_t* q = bam_init1();

    if (align.chromName.empty())
    {
        q->core.tid = -1;
    }
    else
    {
        q->core.tid = bam_name2id(bamHeader_.get(), align.chromName.c_str());
        if (q->core.tid == -1)
        {
            throw std::logic_error("Unknown contig name " + align.chromName);
        }
    }
    q->core.pos = align.pos;
    q->core.mtid = -1;
    q->core.mpos = -1;
    q->core.flag = BAM_FUNMAP;
    if (align.isPaired)
    {
        q->core.flag += BAM_FPAIRED + BAM_FMUNMAP;
        q->core.flag += align.isMate1 ? BAM_FREAD1 : BAM_FREAD2;
    }
    q->core.l_qname = align.fragmentName.length() + 1; // +1 includes the tailing '\0'
    q->core.l_qseq = align.sequence.length();
    q->core.n_cigar = 0; // we have no cigar sequence

    //`q->data` structure: qname-cigar-seq-qual-aux
    int seqQualLength = (int)(1.5 * align.sequence.length() + (align.sequence.length() % 2 != 0));
    q->l_data = q->core.l_qname + seqQualLength;
    q->m_data = q->l_data;
    kroundup32(q->m_data);
    q->data = (uint8_t*)realloc(q->data, q->m_data);
    memcpy(q->data, align.fragmentName.c_str(), q->core.l_qname); // first set qname
    uint8_t* s = bam_get_seq(q);
    for (int i = 0; i < q->core.l_qseq; ++i)
    {
        bam1_seq_seti(s, i, seq_nt16_table[(unsigned char)align.sequence[i]]);
    }
    s = bam_get_qual(q);
    if (!align.BaseQualities.empty() && align.BaseQualities.size() != align.sequence.length())
    {
        throw std::logic_error("Mismatched sequence and quality lengths");
    }
    for (unsigned i = 0; i < align.sequence.length(); ++i)
    {
        s[i] = align.BaseQualities.empty() ? 0xFF : align.BaseQualities[i];
    }
    if (!align.graphCigar.empty())
    {
        string data(align.graphCigar);
        bam_aux_append(
            q, BamWriter::graphCigarBamTag.c_str(), 'Z', align.graphCigar.length() + 1,
            reinterpret_cast<uint8_t*>(&data[0]));
    }
    if (bam_write1(filePtr_->fp.bgzf, q) == 0)
    {
        throw std::logic_error("Cannot write alignment");
    }
    bam_destroy1(q);
}

BamAlignment BamWriter::makeAlignment(
    string const& fragmentName, Sequence const& sequence, std::vector<int> const& qualities,
    BamWriter::PairingInfo pairing, string const& graphAlign) const
{
    BamAlignment align;
    align.fragmentName = fragmentName;
    align.sequence = sequence;
    align.BaseQualities == qualities;
    align.graphCigar = graphAlign;
    align.isPaired = (pairing != BamWriter::PairingInfo::Unpaired);
    align.isMate1 = (pairing == BamWriter::PairingInfo::FirstMate);
    return align;
}

BamAlignment BamWriter::makeAlignment(
    GraphReferenceMapping const& refMap, string const& fragmentName, Sequence const& sequence,
    std::vector<int> const& qualities, BamWriter::PairingInfo pairing, GraphAlignment const& graphAlign) const
{
    BamAlignment bamAlign;
    bamAlign.fragmentName = fragmentName;
    bamAlign.sequence = sequence;
    bamAlign.BaseQualities = qualities;
    bamAlign.isPaired = (pairing != BamWriter::PairingInfo::Unpaired);
    bamAlign.isMate1 = (pairing == BamWriter::PairingInfo::FirstMate);
    auto refPos = refMap.map(graphAlign.path());
    if (refPos)
    {
        bamAlign.chromName = refPos->contig;
        bamAlign.pos = refPos->start;
    }
    std::stringstream graphCigar;
    graphCigar << graphAlign;
    bamAlign.graphCigar = graphCigar.str();
    return bamAlign;
}
}