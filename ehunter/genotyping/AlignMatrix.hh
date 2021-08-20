//
// ExpansionHunter
// Copyright 2016-2020 Illumina, Inc.
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

#include <iosfwd>
#include <limits>
#include <vector>

#include "graphalign/GraphAlignment.hh"

#include "genotyping/StrAlign.hh"

namespace ehunter
{
namespace strgt
{

class AlignMatrix
{
public:
    explicit AlignMatrix(int strNode);
    int numReads() const { return alignScoreMatrix_.size(); }
    void add(const graphtools::GraphAlignment& read, const graphtools::GraphAlignment& mate);
    void remove(int readIndex);
    StrAlign getAlign(int readIndex, int alleleSize) const;
    StrAlign getBestAlign(int readIndex) const;
    const std::vector<std::vector<StrAlign>>& matrix() const { return alignScoreMatrix_; }
    int getMaxMotifCount() const;

    friend void addIrrPairsIfPossibleExpansion(int maxMotifsInRead, AlignMatrix& alignMatrix, int numIrrPairs);

private:
    void add(const graphtools::GraphAlignment& graphAlign);
    int strNode_;
    ConsistentAlignmentCalculator alignmentCalculator_;
    std::vector<StrAlign> bestAlignsByRead_;
    std::vector<std::vector<StrAlign>> alignScoreMatrix_;
};

std::ostream& operator<<(std::ostream& out, const AlignMatrix& matrix);

void addIrrPairsIfPossibleExpansion(int maxMotifsInRead, AlignMatrix& alignMatrix, int numIrrPairs);

}
}
