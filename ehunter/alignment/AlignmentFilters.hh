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

#include <string>

#include <boost/optional.hpp>

#include "graphalign/GraphAlignment.hh"

namespace ehunter
{

/**
 * Checks if a read pair is likely to have originated in the alignment region
 *
 * The check is performed by verifying that the alignment score to non-repeat nodes (combined for both mates) is
 * sufficiently high.
 *
 * @param readAlignment: Alignment of a read
 * @param mateAlignment: Alignment of read's mate
 * @param kMinNonRepeatAlignmentScore: Score threshold
 * @return true if the alignment score to non-repeat nodes exceeds the threshold
 */
bool checkIfLocallyPlacedReadPair(
    boost::optional<graphtools::GraphAlignment> readAlignment,
    boost::optional<graphtools::GraphAlignment> mateAlignment, int kMinNonRepeatAlignmentScore);

// Checks if alignment upstream of a given node is high quality
bool checkIfUpstreamAlignmentIsGood(graphtools::NodeId nodeId, graphtools::GraphAlignment alignment);

// Checks if alignment downstream of a given node is high quality
bool checkIfDownstreamAlignmentIsGood(graphtools::NodeId nodeId, graphtools::GraphAlignment alignment);

bool checkIfPassesAlignmentFilters(const graphtools::GraphAlignment& alignment);

}
