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

#include "core/Reference.hh"

#include <algorithm>
#include <memory>
#include <stdexcept>

using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

FastaReference::FastaReference(const string& referencePath, const ReferenceContigInfo& contigInfo)
    : referencePath_(referencePath)
    , fastaContigInfo_({})
    , bamHeaderContigInfo_(contigInfo)
{
    htsFastaIndexPtr_ = fai_load(referencePath_.c_str());

    std::vector<std::pair<std::string, int64_t>> internalNamesAndSizes;

    for (int contigIndex = 0; contigIndex != faidx_nseq(htsFastaIndexPtr_); ++contigIndex)
    {
        const char* sequenceName = faidx_iseq(htsFastaIndexPtr_, contigIndex);
        int64_t sequenceLength = faidx_seq_len(htsFastaIndexPtr_, sequenceName);
        internalNamesAndSizes.emplace_back(sequenceName, sequenceLength);
    }

    fastaContigInfo_ = ReferenceContigInfo(internalNamesAndSizes);
    assertConsistency();
}

void FastaReference::assertConsistency() const
{
    // This is a stub for a function that should throw an error unless all of the following holds
    // 1. Every reference sequence (used in the catalog) has an ID from the BAM header (up to 'chr' mismatches)
    // 2. Same name refers to same chromosome (coordinate system)
}

FastaReference::~FastaReference() { fai_destroy(htsFastaIndexPtr_); }

string FastaReference::getSequence(const string& contigName, int64_t start, int64_t end) const
{
    const int contigIndex = fastaContigInfo_.getContigId(contigName);
    const char* contigNamePtr = faidx_iseq(htsFastaIndexPtr_, contigIndex);

    int extractedLength;
    // This htslib function is 0-based closed but our coordinates are half open
    char* sequencePtr = faidx_fetch_seq(htsFastaIndexPtr_, contigNamePtr, start, end - 1, &extractedLength);

    if (!sequencePtr || extractedLength < 0 || extractedLength < end - start)
    {
        const string encoding(contigName + ":" + to_string(start) + "-" + to_string(end));
        const string message = "Unable to extract " + encoding + " from " + referencePath_;
        throw std::runtime_error(message);
    }

    string sequence(sequencePtr);
    free(sequencePtr);
    std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);

    return sequence;
}

string FastaReference::getSequence(const GenomicRegion& region) const
{
    return getSequence(bamHeaderContigInfo_.getContigName(region.contigIndex()), region.start(), region.end());
}

}
