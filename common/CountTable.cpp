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

#include "common/CountTable.hh"

using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

int32_t CountTable::countOf(int32_t element) const
{
    if (elements_to_counts_.find(element) == elements_to_counts_.end())
    {
        return 0;
    }
    return elements_to_counts_.at(element);
}

void CountTable::setCountOf(int32_t element, int32_t count)
{
    elements_to_counts_[element] = count;
    if (count == 0)
    {
        elements_to_counts_.erase(element);
    }
}

void CountTable::incrementCountOf(int32_t element) { ++elements_to_counts_[element]; }

vector<int32_t> CountTable::getElementsWithNonzeroCounts() const
{
    vector<int32_t> elements;
    for (const auto& element_count : elements_to_counts_)
    {
        elements.push_back(element_count.first);
    }

    return elements;
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

}
