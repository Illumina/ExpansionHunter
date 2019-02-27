//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
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
//

#pragma once

#include <iostream>
#include <list>
#include <string>

#include "graphalign/GraphAlignment.hh"

namespace ehunter
{

/**
 * Softclips unreliable prefix of an alignment
 *
 * To determine if a prefix of the alignment is unreliable, the prefix is realigned along all valid alternate paths. The
 * alignment is then shrank to the point where high-scoring prefix alignments diverge.
 *
 * @param referenceLength: Reference length of the prefix to evaluate
 * @param query: Query sequence corresponding to the alignment
 * @param alignment: Any graph alignment
 */
void shrinkUncertainPrefix(int referenceLength, const std::string& query, graphtools::GraphAlignment& alignment);

/**
 * Softclips unreliable suffix of an alignment
 *
 * Works identically to shrinkUncertainPrefix but for suffixes
 *
 * @param referenceLength: Reference length of the suffix to evaluate
 * @param query: Query sequence corresponding to the alignment
 * @param alignment: Any graph alignment
 */
void shrinkUncertainSuffix(int referenceLength, const std::string& query, graphtools::GraphAlignment& alignment);

}
