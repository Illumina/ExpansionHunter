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

#include "common/ReferenceContigInfo.hh"

#include <memory>
#include <stdexcept>

using std::pair;
using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{
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
    const auto entry = nameToIndex_.find(contigName);
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