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

#include "common/GenomicRegion.hh"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>

#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using std::istream;
using std::ostream;
using std::string;
using std::vector;

using boost::lexical_cast;
using boost::algorithm::is_any_of;
using boost::algorithm::split;

namespace ehunter
{

Region::Region(const std::string chrom, int64_t start, int64_t end)
    : chrom_(chrom)
    , start_(start)
    , end_(end)
{
}

Region::Region(const std::string encoding)
{
    vector<string> components;
    boost::algorithm::split(components, encoding, is_any_of(":-"));

    if (components.size() != 3)
    {
        throw std::logic_error("Unexpected range format: " + encoding);
    }

    chrom_ = components[0];
    start_ = lexical_cast<int64_t>(components[1]);
    end_ = lexical_cast<int64_t>(components[2]);
}

bool Region::operator<(const Region& other) const
{
    if (chrom_ != other.chrom_)
    {
        return chrom_ < other.chrom_;
    }

    if (start_ != other.start_)
    {
        return start_ < other.start_;
    }

    return end_ < other.end_;
}

bool Region::Overlaps(const Region& other) const
{
    if (chrom_ != other.chrom_)
    {
        return false;
    }

    const int64_t leftBound = start_ > other.start_ ? start_ : other.start_;
    const int64_t rightBound = end_ < other.end_ ? end_ : other.end_;

    return leftBound <= rightBound;
}

int64_t Region::Distance(const Region& other) const
{
    if (chrom_ != other.chrom_)
    {
        return std::numeric_limits<int64_t>::max();
    }

    if (end_ < other.start_)
    {
        return other.start_ - end_;
    }

    if (other.end_ < start_)
    {
        return start_ - other.end_;
    }

    return 0;
}

vector<Region> merge(vector<Region> regions, int maxMergeDist)
{
    if (regions.empty())
    {
        return regions;
    }

    std::sort(regions.begin(), regions.end());

    Region mergedRegion = regions.front();
    vector<Region> mergedRegions;

    for (const auto& currentRegion : regions)
    {
        if (currentRegion.Distance(mergedRegion) <= maxMergeDist)
        {
            const int64_t furthestEnd = std::max<int64_t>(mergedRegion.end(), currentRegion.end());
            mergedRegion.setEnd(furthestEnd);
        }
        else
        {
            mergedRegions.push_back(mergedRegion);
            mergedRegion = currentRegion;
        }
    }

    if (mergedRegions.empty() || (mergedRegions.back() != mergedRegion))
    {
        mergedRegions.push_back(mergedRegion);
    }

    return mergedRegions;
}

const string Region::ToString() const
{
    std::ostringstream out;
    out << *this;
    return out.str();
}

// Returns the range extended by flankSize upstream and downstream.
// NOTE: The right boundary of the extended region may stick past chromosome
// end.
Region Region::extend(int length) const
{
    const int64_t new_start = start_ >= length ? (start_ - length) : 0;
    const int64_t new_end = end_ + length;
    return Region(chrom_, new_start, new_end);
}

std::ostream& operator<<(std::ostream& out, const Region& region)
{
    out << region.chrom_ << ":" << region.start_ << "-" << region.end_;
    return out;
}

}
