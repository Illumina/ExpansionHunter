//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "io/VcfWriterHelpers.hh"

#include <boost/algorithm/string/join.hpp>

#include <sstream>
#include <vector>

using std::deque;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::stringstream;
using std::to_string;
using std::vector;

namespace ehunter
{

template <typename T> static string encodeSampleFields(const deque<T>& fieldRecords)
{
    vector<string> fieldRecordEncodings;

    for (const auto& record : fieldRecords)
    {
        stringstream stringStream;
        stringStream << record;
        fieldRecordEncodings.emplace_back(stringStream.str());
    }

    return boost::algorithm::join(fieldRecordEncodings, "/");
}

string VcfAlleleFields::encode() const
{
    vector<string> encoding;
    encoding.push_back(encodeSampleFields(genotype_));
    encoding.push_back(encodeSampleFields(sources_));
    encoding.push_back(encodeSampleFields(alleleSizes_));
    encoding.push_back(encodeSampleFields(confidenceIntervals_));
    encoding.push_back(encodeSampleFields(spanningReadCounts_));
    encoding.push_back(encodeSampleFields(flankingReadCounts_));
    encoding.push_back(encodeSampleFields(repeatReadCounts_));

    return boost::algorithm::join(encoding, ":");
}

void VcfAlleleFields::addAlleleInfo(
    int alleleSize, ReadType source, NumericInterval confidenceInterval, int spanningReadCount, int flankingReadCount,
    int repeatReadCount)
{
    if (alleleSize == referenceSize_)
    {
        addRefAlleleInfo(alleleSize, source, confidenceInterval, spanningReadCount, flankingReadCount, repeatReadCount);
    }
    else
    {
        addAltAlleleInfo(alleleSize, source, confidenceInterval, spanningReadCount, flankingReadCount, repeatReadCount);
    }
}

void VcfAlleleFields::addRefAlleleInfo(
    int alleleSize, ReadType source, NumericInterval confidenceInterval, int spanningReadCount, int flankingReadCount,
    int repeatReadCount)
{
    genotype_.push_front(0);
    alleleSizes_.push_front(alleleSize);
    sources_.push_front(source);
    confidenceIntervals_.push_front(confidenceInterval);
    spanningReadCounts_.push_front(spanningReadCount);
    flankingReadCounts_.push_front(flankingReadCount);
    repeatReadCounts_.push_front(repeatReadCount);
}

void VcfAlleleFields::addAltAlleleInfo(
    int alleleSize, ReadType source, NumericInterval confidenceInterval, int spanningReadCount, int flankingReadCount,
    int repeatReadCount)
{
    const int previousAlleleSize = alleleSizes_.empty() ? -1 : alleleSizes_.back();
    if (alleleSize < previousAlleleSize)
    {
        throw std::logic_error(
            "Allele of size " + to_string(alleleSize) + " cannot follow allele of size "
            + to_string(previousAlleleSize));
    }

    int haplotypeNum = 1;
    if (!genotype_.empty())
    {
        const int previousAlleleSize = alleleSizes_.back();
        haplotypeNum = previousAlleleSize == alleleSize ? genotype_.back() : genotype_.back() + 1;
    }

    genotype_.push_back(haplotypeNum);
    alleleSizes_.push_back(alleleSize);
    sources_.push_back(source);
    confidenceIntervals_.push_back(confidenceInterval);
    spanningReadCounts_.push_back(spanningReadCount);
    flankingReadCounts_.push_back(flankingReadCount);
    repeatReadCounts_.push_back(repeatReadCount);
}

}
