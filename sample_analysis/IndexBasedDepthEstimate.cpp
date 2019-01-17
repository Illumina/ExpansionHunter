//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "sample_analysis/IndexBasedDepthEstimate.hh"

#include <cassert>
#include <iostream>
#include <unordered_set>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>

extern "C"
{
#include "htslib/hts.h"
#include "htslib/sam.h"
}

#include "common/HtsHelpers.hh"

using std::string;
using std::unordered_set;
using namespace boost::accumulators;

namespace ehunter
{

static bool isAutosome(const string& contigName)
{
    static unordered_set<string> autosomeNames
        = { "chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9",  "chr10", "chr11",
            "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
            "1",     "2",     "3",     "4",     "5",     "6",     "7",     "8",     "9",     "10",    "11",
            "12",    "13",    "14",    "15",    "16",    "17",    "18",    "19",    "20",    "21",    "22" };

    auto contigNameIter = autosomeNames.find(contigName);
    if (contigNameIter != autosomeNames.end())
    {
        return true;
    }

    return false;
}

double estimateDepthFromHtsIndex(const std::string& htsFilePath, int readLength)
{
    htsFile* htsFilePtr = sam_open(htsFilePath.c_str(), "r");
    if (!htsFilePtr)
    {
        throw std::runtime_error("Failed to open HTS file " + htsFilePath);
    }

    bam_hdr_t* htsHeaderPtr = sam_hdr_read(htsFilePtr);
    if (!htsHeaderPtr)
    {
        throw std::runtime_error("Failed to load header of " + htsFilePath);
    }

    hts_idx_t* htsIndexPtr = sam_index_load(htsFilePtr, htsFilePath.c_str());
    if (!htsIndexPtr)
    {
        throw std::runtime_error("Failed to load index of " + htsFilePath);
    }

    const auto contigInfo = htshelpers::decodeContigInfo(htsHeaderPtr);

    accumulator_set<double, features<tag::median>> contigDepths;

    for (int contigIndex = 0; contigIndex != contigInfo.numContigs(); ++contigIndex)
    {
        uint64_t numMappedReads, numUnmappedReads;
        hts_idx_get_stat(htsIndexPtr, contigIndex, &numMappedReads, &numUnmappedReads);
        if (isAutosome(contigInfo.getContigName(contigIndex)))
        {
            const int64_t contigLength = contigInfo.getContigSize(contigIndex);
            const double contigDepth = (readLength * numMappedReads) / static_cast<double>(contigLength);
            contigDepths(contigDepth);
        }
    }

    return median(contigDepths);
}

}
