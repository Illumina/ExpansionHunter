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

#pragma once

#include <iostream>
#include <string>
#include <vector>

class Region
{
public:
    friend std::ostream& operator<<(std::ostream& out, const Region& region);

    Region(const std::string chrom, int64_t start, int64_t end);
    Region(const std::string encoding);

    bool operator<(const Region& other) const;

    bool Overlaps(const Region& other) const;
    int64_t Distance(const Region& other) const;

    const std::string& chrom() const { return chrom_; }
    int64_t start() const { return start_; }
    int64_t end() const { return end_; }
    int64_t length() const { return end_ - start_ + 1; }

    void setChrom(const std::string& chrom) { chrom_ = chrom; }
    void setStart(int64_t start) { start_ = start; }
    void setEnd(int64_t end) { end_ = end; }
    bool operator==(const Region& other) const
    {
        return chrom_ == other.chrom_ && start_ == other.start_ && end_ == other.end_;
    }
    bool operator!=(const Region& other) const { return !(*this == other); }

    Region Extend(int length) const;
    const std::string ToString() const;

private:
    std::string chrom_;
    int64_t start_;
    int64_t end_;
};

std::vector<Region> merge(std::vector<Region> regions, int maxMergeDist = 500);
std::ostream& operator<<(std::ostream& out, const Region& region);
