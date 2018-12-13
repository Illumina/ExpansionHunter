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
