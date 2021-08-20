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

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

#include "core/Common.hh"
#include "core/GenomicRegion.hh"
#include "core/Reference.hh"

namespace ehunter
{

class LocusStats
{
public:
    LocusStats(
        AlleleCount alleleCount = AlleleCount::kOne, int meanReadLen = 0, int medianFragLen = 0, double depth = 0)
        : alleleCount_(alleleCount)
        , meanReadLen_(meanReadLen)
        , medianFragLen_(medianFragLen)
        , depth_(depth)
    {
    }

    AlleleCount alleleCount() const { return alleleCount_; }
    int meanReadLength() const { return meanReadLen_; }
    int medianFragLength() const { return medianFragLen_; }
    double depth() const { return depth_; }
    void setDepth(double depth) { depth_ = depth; }

    bool operator==(const LocusStats& other) const;

private:
    AlleleCount alleleCount_;
    int meanReadLen_;
    int medianFragLen_;
    double depth_;
};

std::ostream& operator<<(std::ostream& out, const LocusStats& stats);

// Computes read and coverage statistics for each locus from reads aligning to the flanks
class LocusStatsCalculator
{
public:
    LocusStatsCalculator(ChromType chromType, const graphtools::Graph& graph);

    void inspect(const graphtools::GraphAlignment& readAlign, const graphtools::GraphAlignment& mateAlign);
    void inspectRead(const graphtools::GraphAlignment& readAlign);

    LocusStats estimate(Sex sampleSex);
    void recordReadLen(const graphtools::GraphAlignment& readAlign);

private:
    using AccumulatorStats
        = boost::accumulators::features<boost::accumulators::tag::count, boost::accumulators::tag::mean>;
    using Accumulator = boost::accumulators::accumulator_set<int, AccumulatorStats>;

    void recordFragLen(const graphtools::GraphAlignment& readAlign, const graphtools::GraphAlignment& mateAlign);

    ChromType chromType_;
    Accumulator readLengthAccumulator_;
    Accumulator fragLengthAccumulator_;
    graphtools::NodeId leftFlankId_;
    graphtools::NodeId rightFlankId_;
    int leftFlankLength_;
    int rightFlankLength_;
};

}
