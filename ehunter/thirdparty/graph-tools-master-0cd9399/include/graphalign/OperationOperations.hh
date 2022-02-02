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

#include "graphalign/Operation.hh"

#include <cstdint>
#include <utility>

namespace graphtools
{

// Returns true if a given operation is consistent with the given query and reference sequences
bool checkConsistency(const Operation& operation, const std::string& reference, const std::string& query);

using OperationPair = std::pair<Operation, Operation>;

/**
 * Splits a given operation by reference length
 *
 * @param operation: Any operation that spans over one base of the reference
 * @param prefix_reference_length: Length of the first piece (prefix)
 * @return A pair of operations produced by the split
 */
OperationPair splitByReferenceLength(const Operation& operation, uint32_t prefix_reference_length);
}
