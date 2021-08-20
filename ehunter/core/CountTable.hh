//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
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

#pragma once

#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

namespace ehunter
{

class CountTable
{
public:
    using const_iterator = std::map<int32_t, int32_t>::const_iterator;
    const_iterator begin() const { return elementsToCounts_.begin(); }
    const_iterator end() const { return elementsToCounts_.end(); }

    CountTable() = default;
    explicit CountTable(std::map<int32_t, int32_t> elementsToCounts)
        : elementsToCounts_(std::move(elementsToCounts)) {};
    CountTable(const CountTable& other) = default;

    void clear() { elementsToCounts_.clear(); }

    int32_t countOf(int32_t element) const;
    void incrementCountOf(int32_t element, int32_t increment = 1);
    void setCountOf(int32_t element, int32_t count);
    std::vector<int32_t> getElementsWithNonzeroCounts() const;

    CountTable& operator=(const CountTable& other);
    bool operator==(const CountTable& other) const { return elementsToCounts_ == other.elementsToCounts_; }

private:
    std::map<int32_t, int32_t> elementsToCounts_;
};

std::ostream& operator<<(std::ostream& out, const CountTable& count_table);

/**
 * Collapses the elements above a specified bound; the counts of collapsed elements are summed together.
 *
 * @param countTable: Any count table
 * @param upperBound: Elements above this bound are collapsed
 * @return Table where counts of all elements exceeding the bound are removed and their counts are added to bound's
 * count
 */
CountTable collapseTopElements(const CountTable& countTable, int upperBound);

}
