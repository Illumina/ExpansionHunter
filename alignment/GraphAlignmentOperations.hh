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

#include "graphalign/GraphAlignment.hh"
#include "graphalign/LinearAlignmentParameters.hh"

namespace ehunter
{

/**
 * Adds softclips to the ends of the alignment
 *
 * @param alignment: any alignment
 * @param leftSoftclipLen: length of the left softclip to add
 * @param rightSoftclipLen: length of the right softclip to add
 * @return extended alignment
 *
 * Example:
 *  GraphAlignment alignment = decodeGraphAlignment(1, "0[3M]1[3M]", &graph);
 *  GraphAlignment extendedAlignment = extendWithSoftclip(alignment, 5, 4);
 *  // extendedAlignment == decodeGraphAlignment(1, "0[5S3M]1[3M4S]", &graph);
 */
graphtools::GraphAlignment
extendWithSoftclip(const graphtools::GraphAlignment& alignment, int leftSoftclipLen, int rightSoftclipLen);

int getNumNonrepeatMatchesUpstream(graphtools::NodeId nodeId, graphtools::GraphAlignment alignment);

int getNumNonrepeatMatchesDownstream(graphtools::NodeId nodeId, graphtools::GraphAlignment alignment);

int scoreAlignmentToNonloopNodes(
    graphtools::GraphAlignment alignment, LinearAlignmentParameters parameters = LinearAlignmentParameters());

int countFullOverlaps(graphtools::NodeId nodeId, graphtools::GraphAlignment alignment);

graphtools::GraphAlignment computeCanonicalAlignment(const std::list<graphtools::GraphAlignment>& alignments);

}
