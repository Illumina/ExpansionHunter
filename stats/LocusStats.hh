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

#include <iostream>
#include <string>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/optional.hpp>

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

#include "common/Common.hh"
#include "common/GenomicRegion.hh"
#include "common/Reference.hh"

namespace ehunter
{

class LocusStats
{
public:
    LocusStats(AlleleCount alleleCount, int meanReadLength, double depth)
        : alleleCount_(alleleCount)
        , meanReadLength_(meanReadLength)
        , depth_(depth)
    {
    }

    AlleleCount alleleCount() const { return alleleCount_; }
    int meanReadLength() const { return meanReadLength_; }
    double depth() const { return depth_; }

    bool operator==(const LocusStats& other) const;

private:
    AlleleCount alleleCount_;
    int meanReadLength_;
    double depth_;
};

std::ostream& operator<<(std::ostream& out, const LocusStats& stats);
}
