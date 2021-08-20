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

#include "io/BamletWriter.hh"

#include <vector>

#include <boost/algorithm/string/join.hpp>

// cppcheck-suppress missingInclude
#include "htslib/bgzf.h"
// cppcheck-suppress missingInclude
#include "htslib/kseq.h"
// cppcheck-suppress missingInclude
#include "htslib/sam.h"

using graphtools::GraphAlignment;
using graphtools::GraphReferenceMapping;
using graphtools::ReferenceInterval;
using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

static GraphReferenceMapping generateMapping(const ReferenceContigInfo& contigInfo, const LocusSpecification& locusSpec)
{
    GraphReferenceMapping mapping(&locusSpec.regionGraph());

    for (const auto& nodeAndRegion : locusSpec.referenceProjectionOfNodes())
    {
        auto nodeId = nodeAndRegion.first;
        const auto& region = nodeAndRegion.second;
        const auto& regionEncoding = encode(contigInfo, region);
        ReferenceInterval referenceInterval = ReferenceInterval::parseRegion(regionEncoding);
        mapping.addMapping(nodeId, referenceInterval);
    }

    return mapping;
}

// Set a base in the bit-backed read sequence representation of a BAM record. (4 bits per base).
// Sets position i to base c in sequence s.
#define bam1_seq_seti(s, i, c) ((s)[(i) >> 1] = ((s)[(i) >> 1] & 0xf << (((i)&1) << 2)) | (c) << ((~(i)&1) << 2))

BamletWriter::BamletWriter(
    const string& bamletPath, const ReferenceContigInfo& contigInfo, const RegionCatalog& regionCatalog)
    : filePtr_(hts_open(bamletPath.c_str(), "wb"), hts_close)
    , bamHeader_(bam_hdr_init(), bam_hdr_destroy)
    , contigInfo_(contigInfo)
    , writeThread_(&BamletWriter::writeHtsAlignments, this)
{
    for (const auto& locusSpec : regionCatalog)
    {
        graphReferenceMappings_.emplace(locusSpec.locusId(), generateMapping(contigInfo, locusSpec));
    }

    writeHeader();
}

BamletWriter::~BamletWriter()
{
    writeQueue_.push(nullptr);
    writeThread_.join();
}

void BamletWriter::writeHeader()
{
    const string initHeader = "@HD\tVN:1.4\tSO:unknown\n";
    bamHeader_->l_text = strlen(initHeader.c_str());
    bamHeader_->text = strdup(initHeader.c_str());
    bamHeader_->n_targets = contigInfo_.numContigs();

    // All this memory gets freed by the header (bam_hdr_destroy) using free
    bamHeader_->target_len = (uint32_t*)calloc(contigInfo_.numContigs(), sizeof(uint32_t));
    bamHeader_->target_name = (char**)calloc(contigInfo_.numContigs(), sizeof(char*));

    for (int index = 0; index != contigInfo_.numContigs(); ++index)
    {
        const string& contigName = contigInfo_.getContigName(index);
        bamHeader_->target_name[index] = (char*)malloc(contigName.length() + 1);
        memcpy(bamHeader_->target_name[index], contigName.c_str(), contigName.length() + 1);
        bamHeader_->target_len[index] = contigInfo_.getContigSize(index);
    }

    if (bam_hdr_write(filePtr_->fp.bgzf, bamHeader_.get()) != 0)
    {
        throw std::logic_error("Failed to write header");
    }
}

void BamletWriter::write(
    const string& locusId, const string& fragmentName, const string& query, bool isFirstMate, bool isReversed,
    bool isMateReversed, const GraphAlignment& alignment)
{
    const GraphReferenceMapping& referenceMapping = graphReferenceMappings_.at(locusId);
    auto optionalReferenceInterval = referenceMapping.map(alignment.path());

    if (optionalReferenceInterval)
    {
        write(*optionalReferenceInterval, fragmentName, query, isFirstMate, isReversed, isMateReversed, alignment);
    }
    else
    {
        ReferenceInterval interval("", -1, -1);
        write(interval, fragmentName, query, isFirstMate, isReversed, isMateReversed, alignment);
    }
}

// TODO: Consider moving to graph tools
static string summarizeAlignment(const GraphAlignment& alignment)
{
    const vector<string> customTagComponents
        = { alignment.path().graphRawPtr()->graphId, to_string(alignment.path().startPosition()),
            alignment.generateCigar() };

    return boost::algorithm::join(customTagComponents, ",");
}

static const string graphAlignmentBamTag = "XG";

static vector<int> extractQualityScores(const string& query)
{
    const int kLowQualityScore = 0;
    const int kHighQualityScore = 40;

    vector<int> qualities;
    qualities.resize(query.length());

    for (int index = 0; index != static_cast<int>(query.size()); ++index)
    {
        qualities[index] = isupper(query[index]) ? kHighQualityScore : kLowQualityScore;
    }

    return qualities;
}

void BamletWriter::write(
    const ReferenceInterval& interval, const string& fragmentName, const string& query, bool isFirstMate,
    bool isReversed, bool isMateReversed, const GraphAlignment& alignment)
{

    bam1_t* htsAlignmentPtr = bam_init1();

    if (interval.contig.empty())
    {
        htsAlignmentPtr->core.tid = -1;
    }
    else
    {
        htsAlignmentPtr->core.tid = bam_name2id(bamHeader_.get(), interval.contig.c_str());
        if (htsAlignmentPtr->core.tid == -1)
        {
            throw std::logic_error("Unknown contig name " + interval.contig);
        }
    }

    htsAlignmentPtr->core.pos = interval.start;
    htsAlignmentPtr->core.mtid = -1;
    htsAlignmentPtr->core.mpos = -1;
    htsAlignmentPtr->core.flag = BAM_FUNMAP;

    htsAlignmentPtr->core.flag += BAM_FPAIRED + BAM_FMUNMAP;

    if (isReversed)
        htsAlignmentPtr->core.flag += BAM_FREVERSE;
    if (isMateReversed)
        htsAlignmentPtr->core.flag += BAM_FMREVERSE;

    htsAlignmentPtr->core.flag += isFirstMate ? BAM_FREAD1 : BAM_FREAD2;

    const vector<int> qualities = extractQualityScores(query);

    htsAlignmentPtr->core.l_qname = fragmentName.length() + 1; // +1 includes the tailing '\0'
    htsAlignmentPtr->core.l_qseq = query.length();
    htsAlignmentPtr->core.n_cigar = 0; // we have no cigar sequence

    //`q->data` structure: qname-cigar-seq-qual-aux
    int seqQualLength = (int)(1.5 * query.length() + (query.length() % 2 != 0));
    htsAlignmentPtr->l_data = htsAlignmentPtr->core.l_qname + seqQualLength;
    htsAlignmentPtr->m_data = htsAlignmentPtr->l_data;
    kroundup32(htsAlignmentPtr->m_data);
    htsAlignmentPtr->data = (uint8_t*)realloc(htsAlignmentPtr->data, htsAlignmentPtr->m_data);
    memcpy(htsAlignmentPtr->data, fragmentName.c_str(), htsAlignmentPtr->core.l_qname); // first set qname

    uint8_t* htsSequencePtr = bam_get_seq(htsAlignmentPtr);
    for (int index = 0; index < htsAlignmentPtr->core.l_qseq; ++index)
    {
        bam1_seq_seti(htsSequencePtr, index, seq_nt16_table[(unsigned char)query[index]]);
    }

    htsSequencePtr = bam_get_qual(htsAlignmentPtr);
    if (!qualities.empty() && qualities.size() != query.length())
    {
        throw std::logic_error("Mismatched sequence and quality lengths");
    }

    for (unsigned index = 0; index < query.length(); ++index)
    {
        htsSequencePtr[index] = qualities.empty() ? 0xFF : qualities[index];
    }

    string alignmentEncoding = summarizeAlignment(alignment);
    bam_aux_append(
        htsAlignmentPtr, graphAlignmentBamTag.c_str(), 'Z', alignmentEncoding.length() + 1,
        reinterpret_cast<uint8_t*>(&alignmentEncoding[0]));

    writeQueue_.push(htsAlignmentPtr);
}

void BamletWriter::writeHtsAlignments()
{
    bam1_t* htsAlignmentPtr(nullptr);
    while (true)
    {
        writeQueue_.pop(htsAlignmentPtr);
        if (not htsAlignmentPtr)
            return;
        if (bam_write1(filePtr_->fp.bgzf, htsAlignmentPtr) == 0)
        {
            throw std::logic_error("Cannot write alignment");
        }
        bam_destroy1(htsAlignmentPtr);
    }
}

}
