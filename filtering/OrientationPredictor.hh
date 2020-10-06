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

#include <cstdint>
#include <iostream>
#include <memory>

#include "graphalign/GraphAlignment.hh"
#include "graphalign/KmerIndex.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphOperations.hh"

namespace ehunter
{

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
    OrientationPredictor(const graphtools::Graph* graphRawPtr)
        : kmerLength_(10)
        , minKmerMatchesToPass_(3)
        , kmerIndex_(*graphRawPtr, kmerLength_)
        , reverseComplementedGraph_(graphtools::reverseGraph(*graphRawPtr, true))
        , kmerIndexForReverseComplementedGraph_(reverseComplementedGraph_, kmerLength_)
    {
    }

    OrientationPrediction predict(const std::string& query) const;

private:
    int32_t kmerLength_;
    int32_t minKmerMatchesToPass_;
    graphtools::KmerIndex kmerIndex_;
    const graphtools::Graph reverseComplementedGraph_;
    graphtools::KmerIndex kmerIndexForReverseComplementedGraph_;
};

}
