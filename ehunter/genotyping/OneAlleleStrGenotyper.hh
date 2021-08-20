//
// ExpansionHunter
// Copyright 2016-2020 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Konrad Scheffler <kscheffler@illumina.com>,
//         Felix Schlesinger <fschlesinger@illumina.com>
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

#include <unordered_set>

#include "genotyping/AlignMatrix.hh"
#include "genotyping/FragLogliks.hh"
#include "genotyping/RepeatGenotype.hh"

namespace ehunter
{
namespace strgt
{

class OneAlleleGenotyper
{
public:
    OneAlleleGenotyper(int motifLen, std::vector<double> topFragLogliks, FragLogliks* fragLogliksPtr)
        : motifLen_(motifLen)
        , topFragLogliks_(std::move(topFragLogliks))
        , fragLogliks_(*fragLogliksPtr)
    {
    }

    RepeatGenotype genotype(const std::unordered_set<int>& alleleSizeCandidates);

private:
    RepeatGenotype getMostLikelyGenotype(const std::unordered_set<int>& alleleSizeCandidates);
    double getAlleleLoglik(int motifCount);

    int motifLen_;
    std::vector<double> topFragLogliks_;
    FragLogliks& fragLogliks_;
};

}
}