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
