//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
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

#include "graphalign/OperationOperations.hh"

#include <sstream>
#include <string>

#include "graphutils/BaseMatching.hh"

using std::string;

namespace graphtools
{
bool checkConsistency(const Operation& operation, const std::string& reference, const std::string& query)
{
    const bool is_query_full_length = query.length() == operation.length();
    const bool is_ref_full_length = reference.length() == operation.length();

    switch (operation.type())
    {
    case OperationType::kMatch:
        if (is_query_full_length && checkIfReferenceAndQuerySequencesMatch(reference, query))
            return true;
        return false;

    case OperationType::kMismatch:
        if (is_query_full_length && query.length() == reference.length())
        {
            bool found_matching_base = false;
            for (size_t index = 0; index != query.length(); ++index)
            {
                if (checkIfReferenceBaseMatchesQueryBase(reference[index], query[index]))
                {
                    found_matching_base = true;
                }
            }
            if (!found_matching_base)
            {
                return true;
            }
        }
        return false;

    case OperationType::kMissingBases:
        if (query.length() == reference.length() && is_query_full_length)
        {
            bool found_non_n_base_in_query = false;
            for (size_t index = 0; index != query.length(); ++index)
            {
                if (query[index] != 'N')
                {
                    found_non_n_base_in_query = true;
                }
            }
            if (!found_non_n_base_in_query)
            {
                return true;
            }
        }
        return false;

    case OperationType::kDeletionFromRef:
        if (query.empty() && !reference.empty() && is_ref_full_length)
            return true;
        return false;

    case OperationType::kInsertionToRef:
        if (!query.empty() && reference.empty() && is_query_full_length)
            return true;
        return false;

    case OperationType::kSoftclip:
        if (!query.empty() && reference.empty() && is_query_full_length)
            return true;
        return false;
    }

    return false;
}

OperationPair splitByReferenceLength(const Operation& operation, uint32_t prefix_reference_length)
{
    if (prefix_reference_length == 0 || operation.referenceLength() <= prefix_reference_length)
    {
        std::ostringstream os;
        os << operation;
        const string msg = os.str() + " cannot be split by reference length " + std::to_string(prefix_reference_length);
        throw std::logic_error(msg);
    }

    uint32_t suffix_reference_length = operation.referenceLength() - prefix_reference_length;
    Operation prefix_operation(operation.type(), prefix_reference_length);
    Operation suffix_operation(operation.type(), suffix_reference_length);

    return std::make_pair(prefix_operation, suffix_operation);
}
}
