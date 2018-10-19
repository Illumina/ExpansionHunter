// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
