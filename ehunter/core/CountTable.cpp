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

#include "core/CountTable.hh"

#include <stdexcept>

using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

int32_t CountTable::countOf(int32_t element) const
{
    if (elementsToCounts_.find(element) == elementsToCounts_.end())
    {
        return 0;
    }
    return elementsToCounts_.at(element);
}

void CountTable::setCountOf(int32_t element, int32_t count)
{
    elementsToCounts_[element] = count;
    if (count == 0)
    {
        elementsToCounts_.erase(element);
    }
}

void CountTable::incrementCountOf(int32_t element, int32_t increment)
{
    if (increment <= 0)
    {
        throw std::logic_error("CountTables require positive increments");
    }

    elementsToCounts_[element] += increment;
}

vector<int32_t> CountTable::getElementsWithNonzeroCounts() const
{
    vector<int32_t> elements;
    for (const auto& element_count : elementsToCounts_)
    {
        elements.push_back(element_count.first);
    }

    return elements;
}

CountTable& CountTable::operator=(const CountTable& other)
{
    if (this != &other)
    {
        elementsToCounts_ = other.elementsToCounts_;
    }

    return *this;
}

std::ostream& operator<<(std::ostream& out, const CountTable& count_table)
{
    string encoding;

    for (int32_t element : count_table.getElementsWithNonzeroCounts())
    {

        if (!encoding.empty())
        {
            encoding += ", ";
        }

        encoding += "(" + to_string(element) + ", " + to_string(count_table.countOf(element)) + ")";
    }

    if (encoding.empty())
    {
        encoding = "()";
    }

    out << encoding;

    return out;
}

CountTable collapseTopElements(const CountTable& countTable, int upperBound)
{
    if (upperBound < 0)
    {
        throw std::logic_error("CountTables cannot be truncated to negative values");
    }

    CountTable truncatedTable;

    for (auto element : countTable.getElementsWithNonzeroCounts())
    {
        auto count = countTable.countOf(element);

        if (element < upperBound)
        {
            truncatedTable.setCountOf(element, count);
        }
        else
        {
            truncatedTable.incrementCountOf(upperBound, count);
        }
    }

    return truncatedTable;
}

}
