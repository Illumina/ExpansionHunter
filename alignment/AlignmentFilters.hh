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

#include <list>
#include <string>

#include "graphalign/GraphAlignment.hh"

namespace ehunter
{

/**
 * Checks if a read pair originated in the locus defined by the graph
 *
 * Verifies that there is a pair of read/mate alignments with a sufficiently high combined score to non-repeat nodes
 *
 * @param readAlignments: Alignments of a read
 * @param mateAlignments: Alignments of read's mate
 * @param kMinNonRepeatAlignmentScore: Positive score threshold
 * @return true if read/mate alignment with above properties exists
 */
bool checkIfComesFromGraphLocus(
    const std::list<graphtools::GraphAlignment>& readAlignments,
    const std::list<graphtools::GraphAlignment>& mateAlignments, int kMinNonRepeatAlignmentScore);

// Checks if alignment upstream of a given node is high quality
bool checkIfUpstreamAlignmentIsGood(graphtools::NodeId nodeId, const graphtools::GraphAlignment& alignment);

// Checks if alignment downstream of a given node is high quality
bool checkIfDownstreamAlignmentIsGood(graphtools::NodeId nodeId, const graphtools::GraphAlignment& alignment);

bool checkIfPassesAlignmentFilters(const graphtools::GraphAlignment& alignment);

}
