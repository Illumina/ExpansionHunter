//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
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

#include <cstdint>
#include <iostream>
#include <memory>

#include "graphalign/GraphAlignment.hh"
#include "graphalign/KmerIndex.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphOperations.hh"

enum class OrientationPrediction
{
    kAlignsInOriginalOrientation,
    kAlignsInReverseComplementOrientation,
    kDoesNotAlign
};

std::ostream& operator<<(std::ostream& out, OrientationPrediction orientationPrediction);

class OrientationPredictor
{
public:
    OrientationPredictor(int readLength, const graphtools::Graph* graphRawPtr)
        : kmerLength_(pickKmerLength(readLength))
        , minKmerMatchesToPass_(pickKmerMatchesToPass(readLength))
        , kmerIndex_(*graphRawPtr, kmerLength_)
        , reverseComplementedGraph_(graphtools::reverseGraph(*graphRawPtr, true))
        , kmerIndexForReverseComplementedGraph_(reverseComplementedGraph_, kmerLength_)
    {
    }

    static int pickKmerLength(int readLength)
    {
        int kmerLength = readLength / 10.5;
        kmerLength = std::max(kmerLength, 3);

        return kmerLength;
    }

    static int pickKmerMatchesToPass(int readLength)
    {
        if (readLength >= 100)
        {
            return 3;
        }

        return 3;
    }

    OrientationPrediction predict(const std::string& query) const;

private:
    int32_t kmerLength_;
    int32_t minKmerMatchesToPass_;
    graphtools::KmerIndex kmerIndex_;
    const graphtools::Graph reverseComplementedGraph_;
    graphtools::KmerIndex kmerIndexForReverseComplementedGraph_;
};
