//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include <iostream>
#include <list>
#include <string>

#include "graphalign/GraphAlignment.hh"

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
