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

#include "core/ReferenceContigInfo.hh"

#include <memory>
#include <stdexcept>

using std::pair;
using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

namespace
{
// Removes "chr" prefix from contig names that contain it; adds it to contigs that don't
string generateAlternativeContigName(const string& originalName)
{
    if (originalName.length() > 3 && originalName.substr(0, 3) == "chr")
    {
        return originalName.substr(3);
    }
    else
    {
        return "chr" + originalName;
    }
}
}

ReferenceContigInfo::ReferenceContigInfo(vector<pair<string, int64_t>> namesAndSizes)
    : namesAndSizes_(std::move(namesAndSizes))
{
    for (int index = 0; index != static_cast<int>(namesAndSizes_.size()); ++index)
    {
        const auto& contigName = namesAndSizes_[index].first;
        nameToIndex_.emplace(std::make_pair(contigName, index));
    }
}

const std::string& ReferenceContigInfo::getContigName(int32_t contigIndex) const
{
    assertValidIndex(contigIndex);
    return namesAndSizes_[contigIndex].first;
}

int64_t ReferenceContigInfo::getContigSize(int32_t contigIndex) const
{
    assertValidIndex(contigIndex);
    return namesAndSizes_[contigIndex].second;
}

int32_t ReferenceContigInfo::getContigId(const std::string& contigName) const
{
    auto entry = nameToIndex_.find(contigName);
    if (entry == nameToIndex_.end())
    {
        entry = nameToIndex_.find(generateAlternativeContigName(contigName));
    }

    if (entry == nameToIndex_.end())
    {
        throw std::logic_error("Invalid contig name " + contigName);
    }

    return entry->second;
}

void ReferenceContigInfo::assertValidIndex(int32_t contigIndex) const
{
    if (contigIndex >= static_cast<int32_t>(namesAndSizes_.size()))
    {
        throw std::logic_error("Invalid contig index " + to_string(contigIndex));
    }
}

std::ostream& operator<<(std::ostream& out, const ReferenceContigInfo& contigInfo)
{
    for (int contigIndex = 0; contigIndex != contigInfo.numContigs(); ++contigIndex)
    {
        const auto& contigName = contigInfo.getContigName(contigIndex);
        out << contigName << " -> " << contigIndex << std::endl;
    }

    return out;
}

}
