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

#include <list>
#include <string>

#include "graphalign/GappedAligner.hh"
#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

class GraphAlignmentHeuristicsParameters
{
public:
    GraphAlignmentHeuristicsParameters(int kmerLenForAlignment = 14, int paddingLength = 10, int seedAffixTrimLen = 5)
        : kmerLenForAlignment_(kmerLenForAlignment)
        , paddingLength_(paddingLength)
        , seedAffixTrimLen_(seedAffixTrimLen)
    {
    }
    int kmerLenForAlignment() const { return kmerLenForAlignment_; }
    int paddingLength() const { return paddingLength_; }
    int seedAffixTrimLen() const { return seedAffixTrimLen_; }

private:
    int kmerLenForAlignment_;
    int paddingLength_;
    int seedAffixTrimLen_;
};

class SoftclippingAligner
{
public:
    SoftclippingAligner(
        graphtools::Graph* graphPtr, const std::string& alignerName,
        const GraphAlignmentHeuristicsParameters& alignmentHeuristicsParameters);
    std::list<graphtools::GraphAlignment> align(const std::string& query) const;

private:
    graphtools::GappedGraphAligner aligner_;
};
