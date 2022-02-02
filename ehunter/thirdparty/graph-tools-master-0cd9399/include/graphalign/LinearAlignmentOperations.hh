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

#pragma once

#include <cstdint>
#include <list>
#include <string>
#include <utility>

#include "graphalign/LinearAlignment.hh"
#include "graphutils/SequenceOperations.hh"

namespace graphtools
{

// Checks if a given linear alignment is consistent with the given query and reference sequences
bool checkConsistency(const Alignment& alignment, const std::string& reference, const std::string& query);

/**
 * Splits query and reference into pieces corresponding to operations that the alignment is made of
 *
 * @param alignment: Any linear alignment
 * @param query: Any query sequence
 * @param reference: Any reference sequence
 * @return: List of pairs of sequences corresponding to each operation in the alignment
 */
std::list<StringPair>
getSequencesForEachOperation(const Alignment& alignment, const std::string& reference, const std::string& query);

/**
 * Checks if two alignments are bookended
 *
 * Two alignments are considered bookended if
 *  - Positions of first alignment end and second alignment start are adjacent and
 *  - First alignment does not end in softclipped bases (unless all of its bases are softclipped)
 *  - Second alignment does not start in softclipped bases (unless all of its bases are softclipped)
 *
 * @param first_alignment: Any linear alignment
 * @param second_alignment: Any linear alignment
 * @return true if alignment are bookended
 */
bool checkIfBookended(const Alignment& first_alignment, const Alignment& second_alignment);

/**
 * Merges two bookended alignments into a longer alignment
 *
 * @param first_alignment: Any linear alignment
 * @param second_alignment: Any linear alignment that is bookeneded with the first
 * @return Merged alignment
 */
Alignment mergeAlignments(const Alignment& first_alignment, const Alignment& second_alignment);

/**
 * Calculates alignment score
 *
 * @param alignment: Any linear alignment
 * @param match_score: Score of a match
 * @param mismatch_score: Score of a mismatch
 * @param gap_score: Score of a gap
 * @return Alignment score
 */
int32_t scoreAlignment(const Alignment& alignment, int32_t match_score, int32_t mismatch_score, int32_t gap_score);

// Encodes alignment as a three-row strings where the top corresponds to the query sequence, the bottom to the reference
// sequence, and the middle contains a "|" for each pair of matching bases; gaps are indicated by "-"
std::string prettyPrint(const Alignment& alignment, const std::string& reference, const std::string& query);
}
