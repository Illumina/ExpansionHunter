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
