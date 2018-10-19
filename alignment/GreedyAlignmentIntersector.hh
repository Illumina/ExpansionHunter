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

#include <boost/optional.hpp>

#include "graphalign/GraphAlignment.hh"

class GreedyAlignmentIntersector
{
public:
    GreedyAlignmentIntersector(
        const graphtools::GraphAlignment& firstAlignment, const graphtools::GraphAlignment& secondAlignment)
        : firstAlignment_(firstAlignment)
        , secondAlignment_(secondAlignment)
        , firstPath_(firstAlignment.path())
        , secondPath_(secondAlignment.path())
    {
        initialize();
    }

    boost::optional<graphtools::GraphAlignment> intersect();

private:
    void initialize();
    bool tryAdvancingIndexesToCommonNode();
    bool checkIfAlignmentEndReached(int firstPathIndex, int secondPathIndex);
    bool checkIfCommonNodeIsLoop() const;
    void advanceIndexesToMatchRemainingIterations();
    void advanceIndexesToLastCommonNode();
    void computeIntersectionEnds();
    bool checkIfIntersectionIsConsistent() const;
    boost::optional<graphtools::GraphAlignment> softclipFirstAlignmentToIntersection() const;

    const graphtools::GraphAlignment& firstAlignment_;
    const graphtools::GraphAlignment& secondAlignment_;
    const graphtools::Path& firstPath_;
    const graphtools::Path& secondPath_;

    int nodeIndexOfIntersectionStartOnFirstPath_;
    int nodeIndexOfIntersectionStartOnSecondPath_;
    int nodeIndexOfIntersectionEndOnFirstPath_ = -1;
    int nodeIndexOfIntersectionEndOnSecondPath_ = -1;
    int intersectionStart_ = -1;
    int intersectionEnd_ = -1;
};
