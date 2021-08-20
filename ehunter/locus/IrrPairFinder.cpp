//
// ExpansionHunter
// Copyright 2016-2021 Illumina, Inc.
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

#include "locus/IrrPairFinder.hh"

#include <memory>

using std::string;

namespace ehunter
{
namespace locus
{

IrrPairFinder::IrrPairFinder(string motif)
    : targetMotif_(std::move(motif))
    , purityCalculator_(targetMotif_)
{
}

bool IrrPairFinder::check(const string& read, const string& mate) const
{
    const bool isFirstReadIrr = purityCalculator_.score(read) >= purityCutoff_;
    const bool isSecondReadIrr = purityCalculator_.score(mate) >= purityCutoff_;
    return isFirstReadIrr && isSecondReadIrr;
}

}
}