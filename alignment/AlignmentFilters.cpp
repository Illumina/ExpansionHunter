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

#include "alignment/AlignmentFilters.hh"

#include <list>
#include <vector>

#include "graphalign/GaplessAligner.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphalign/LinearAlignmentOperations.hh"
#include "graphcore/PathOperations.hh"

#include "alignment/GraphAlignmentOperations.hh"

using graphtools::GraphAlignment;
using graphtools::NodeId;
using graphtools::Path;
using std::list;
using std::string;
using std::vector;

namespace ehunter
{

bool checkIfLocallyPlacedReadPair(
    boost::optional<GraphAlignment> readAlignment, boost::optional<GraphAlignment> mateAlignment,
    int kMinNonRepeatAlignmentScore)
{
    int nonRepeatAlignmentScore = 0;

    if (readAlignment)
    {
        nonRepeatAlignmentScore += scoreAlignmentToNonloopNodes(*readAlignment);
    }

    if (mateAlignment)
    {
        nonRepeatAlignmentScore += scoreAlignmentToNonloopNodes(*mateAlignment);
    }

    if (nonRepeatAlignmentScore < kMinNonRepeatAlignmentScore)
    {
        return false;
    }

    return true;
}

bool checkIfUpstreamAlignmentIsGood(NodeId nodeId, GraphAlignment alignment)
{
    const list<int> repeatNodeIndexes = alignment.getIndexesOfNode(nodeId);

    if (repeatNodeIndexes.empty())
    {
        return false;
    }

    const int firstRepeatNodeIndex = repeatNodeIndexes.front();
    int score = 0;
    LinearAlignmentParameters parameters;
    for (int nodeIndex = 0; nodeIndex != firstRepeatNodeIndex; ++nodeIndex)
    {
        score += scoreAlignment(
            alignment[nodeIndex], parameters.matchScore, parameters.mismatchScore, parameters.gapOpenScore);
    }

    const int kScoreCutoff = parameters.matchScore * 8;

    return score >= kScoreCutoff;
}

bool checkIfDownstreamAlignmentIsGood(NodeId nodeId, GraphAlignment alignment)
{
    const list<int> repeatNodeIndexes = alignment.getIndexesOfNode(nodeId);

    if (repeatNodeIndexes.empty())
    {
        return false;
    }

    const int lastRepeatNodeIndex = repeatNodeIndexes.back();
    int score = 0;
    LinearAlignmentParameters parameters;
    for (int nodeIndex = lastRepeatNodeIndex + 1; nodeIndex != static_cast<int>(alignment.size()); ++nodeIndex)
    {
        score += scoreAlignment(
            alignment[nodeIndex], parameters.matchScore, parameters.mismatchScore, parameters.gapOpenScore);
    }

    const int kScoreCutoff = parameters.matchScore * 8;

    return score >= kScoreCutoff;
}

}
