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

#include "graphalign/GappedAligner.hh"
#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

namespace ehunter
{

class SoftclippingAligner
{
public:
    SoftclippingAligner(
        const graphtools::Graph* graphPtr, int kmerLenForAlignment, int paddingLength, int seedAffixTrimLength);
    std::list<graphtools::GraphAlignment>
    align(const std::string& query, graphtools::AlignerSelector& alignerSelector) const;

private:
    graphtools::GappedGraphAligner aligner_;
};

}
