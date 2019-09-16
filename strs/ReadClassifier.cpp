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

#include "strs/ReadClassifier.hh"

using std::vector;

namespace ehunter
{

ReadClassifier::ReadClassifier(vector<GenomicRegion> targetRegions)
    : targetRegions_(std::move(targetRegions))
{
}

ReadClassifier::ReadType ReadClassifier::classify(const MappedRead& read) const
{
    ReadType classification = ReadType::kOfftarget;
    const int64_t readEnd = read.approximateEnd();
    for (const auto& region : targetRegions_)
    {
        if (read.contigIndex() != region.contigIndex())
        {
            continue;
        }

        if (region.start() <= read.pos() && readEnd <= region.end())
        {
            return ReadType::kTarget;
        }

        if (region.start() - kMinOfftargetDistance_ <= read.pos() && readEnd <= region.end() + kMinOfftargetDistance_)
        {
            classification = ReadType::kOther;
        }
    }

    return classification;
}

PairType ReadClassifier::classify(const MappedRead& read, const MappedRead& mate) const
{
    ReadType readType = classify(read);
    ReadType mateType = classify(mate);

    if (readType == ReadType::kTarget || mateType == ReadType::kTarget)
    {
        return PairType::kTarget;
    }

    if (readType == ReadType::kOther || mateType == ReadType::kOther)
    {
        return PairType::kOther;
    }

    return PairType::kOfftarget;
}

std::ostream& operator<<(std::ostream& out, PairType type)
{
    switch (type)
    {
    case PairType::kOther:
        out << "PairType::kOther";
        break;
    case PairType::kTarget:
        out << "PairType::kTarget";
        break;
    case PairType::kOfftarget:
        out << "PairType::kOfftarget";
        break;
    }

    return out;
}

}