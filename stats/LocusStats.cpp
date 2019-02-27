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

#include "stats/LocusStats.hh"

#include <sstream>
#include <stdexcept>

namespace ehunter
{

using boost::optional;

bool LocusStats::operator==(const LocusStats& other) const { return medianReadLength_ == other.medianReadLength_; }

std::ostream& operator<<(std::ostream& out, const LocusStats& stats)
{
    out << "LocusStats(medianReadLength=" << stats.medianReadLength() << ", depth=" << stats.depth() << ")";
    return out;
}

LocusStatsCalculator::LocusStatsCalculator(const graphtools::Graph& graph)
    : readCount_(0)
{
    // As elsewhere in the program, assuming that the fist and last node are flanks
    leftFlankId_ = 0;
    rightFlankId_ = graph.numNodes() - 1;

    leftFlankLength_ = graph.nodeSeq(leftFlankId_).length();
    rightFlankLength_ = graph.nodeSeq(rightFlankId_).length();
}

void LocusStatsCalculator::inspect(const graphtools::GraphAlignment& alignment)
{
    const graphtools::NodeId firstNode = alignment.path().getNodeIdByIndex(0);
    if (firstNode == leftFlankId_ || firstNode == rightFlankId_)
    {
        readLengthAccumulator_(alignment.queryLength());
        ++readCount_;
    }
}

optional<LocusStats> LocusStatsCalculator::estimate() const
{
    if (boost::accumulators::count(readLengthAccumulator_) == 0 || readCount_ == 0)
    {
        return optional<LocusStats>();
    }

    const int meanReadLength = boost::accumulators::mean(readLengthAccumulator_);
    const int numberOfStartPositions = leftFlankLength_ + rightFlankLength_ - meanReadLength;
    const double depth = meanReadLength * (static_cast<double>(readCount_) / numberOfStartPositions);
    const int medianReadLength = boost::accumulators::median(readLengthAccumulator_);
    return LocusStats(medianReadLength, depth);
}

}
