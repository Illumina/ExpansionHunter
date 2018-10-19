//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
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

#include <cstdint>
#include <iostream>
#include <map>
#include <vector>

class CountTable
{
public:
    using const_iterator = std::map<int32_t, int32_t>::const_iterator;
    const_iterator begin() const { return elements_to_counts_.begin(); }
    const_iterator end() const { return elements_to_counts_.end(); }

    CountTable(){};
    explicit CountTable(const std::map<int32_t, int32_t>& elements_to_counts)
        : elements_to_counts_(elements_to_counts){};

    void clear() { elements_to_counts_.clear(); }

    int32_t countOf(int32_t element) const;
    void incrementCountOf(int32_t element);
    void setCountOf(int32_t element, int32_t count);
    std::vector<int32_t> getElementsWithNonzeroCounts() const;

    bool operator==(const CountTable& other) const { return elements_to_counts_ == other.elements_to_counts_; }

private:
    std::map<int32_t, int32_t> elements_to_counts_;
};

std::ostream& operator<<(std::ostream& out, const CountTable& count_table);