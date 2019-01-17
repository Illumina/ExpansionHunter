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

#include "common/Reference.hh"

#include <algorithm>
#include <memory>
#include <stdexcept>

using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

FastaReference::FastaReference(
        const string& referencePath,
        const ReferenceContigInfo& contigInfo)
    : referencePath_(referencePath)
    , contigInfo_(contigInfo)
{
    htsFastaIndexPtr_ = fai_load(referencePath_.c_str());
}

FastaReference::~FastaReference() { fai_destroy(htsFastaIndexPtr_); }

string FastaReference::getSequence(const string& contigName, int64_t start, int64_t end) const
{
    int extractedLength;
    // This htslib function is 0-based closed but our coordinates are half open
    char* sequencePtr = faidx_fetch_seq(htsFastaIndexPtr_, contigName.c_str(), start, end - 1, &extractedLength);

    if (!sequencePtr || extractedLength < 0 || extractedLength < end - start)
    {
        const string encoding(contigName + ":" + to_string(start) + "-" + to_string(end));
        const string message = "Cannot extract " + encoding + " from " + referencePath_
            + "; chromosome names must match exactly (e.g. chr1 and 1 are distinct names) "
            + "and coordinates cannot be past the end of the chromosome";
        throw std::runtime_error(message);
    }

    string sequence("N", extractedLength);
    std::transform(sequencePtr, sequencePtr + extractedLength, sequence.begin(), ::toupper);
    free(sequencePtr);

    return sequence;
}

string FastaReference::getSequence(const GenomicRegion& region) const
{
    return getSequence(contigInfo_.getContigName(region.contigIndex()), region.start(), region.end());
}

}
