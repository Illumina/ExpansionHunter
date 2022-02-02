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

#include "core/LocusStats.hh"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

using graphtools::GraphAlignment;

namespace ehunter
{

bool LocusStats::operator==(const LocusStats& other) const
{
    return alleleCount_ == other.alleleCount_ && meanReadLen_ == other.meanReadLen_
        && medianFragLen_ == other.medianFragLen_ && depth_ == other.depth_;
}

std::ostream& operator<<(std::ostream& out, const LocusStats& stats)
{
    out << "LocusStats(meanReadLength=" << stats.meanReadLength() << ", depth=" << stats.depth() << ")";
    return out;
}

LocusStatsCalculator::LocusStatsCalculator(ChromType chromType, const graphtools::Graph& graph)
    : chromType_(chromType)
{
    // As elsewhere in the program, assuming that the fist and last node are flanks
    leftFlankId_ = 0;
    rightFlankId_ = graph.numNodes() - 1;

    leftFlankLength_ = graph.nodeSeq(leftFlankId_).length();
    rightFlankLength_ = graph.nodeSeq(rightFlankId_).length();
}

void LocusStatsCalculator::inspect(const GraphAlignment& readAlign, const GraphAlignment& mateAlign)
{
    recordReadLen(readAlign);
    recordReadLen(mateAlign);
    recordFragLen(readAlign, mateAlign);
}

void LocusStatsCalculator::inspectRead(const GraphAlignment& readAlign) { recordReadLen(readAlign); }

static AlleleCount determineExpectedAlleleCount(ChromType chromType, Sex sex)
{
    switch (chromType)
    {
    case ChromType::kY:
        return AlleleCount::kOne; // Assume that chrY always has copy number one
    case ChromType::kX:
        return (sex == Sex::kFemale ? AlleleCount::kTwo : AlleleCount::kOne);
    case ChromType::kAutosome:
        return AlleleCount::kTwo;
    }

    return AlleleCount::kTwo; // To remove spurious control reaches end of non-void function warning
}

LocusStats LocusStatsCalculator::estimate(Sex sampleSex)
{
    const int readCount = boost::accumulators::count(readLengthAccumulator_);
    AlleleCount alleleCount = determineExpectedAlleleCount(chromType_, sampleSex);

    if (readCount == 0)
    {
        return { alleleCount, 0, 0, 0.0 };
    }

    const int meanReadLength = boost::accumulators::mean(readLengthAccumulator_);
    const int numberOfStartPositions = leftFlankLength_ + rightFlankLength_ - meanReadLength;
    const double depth = meanReadLength * (static_cast<double>(readCount) / numberOfStartPositions);

    int meanFragLen = 0;
    const int fragCount = boost::accumulators::count(fragLengthAccumulator_);
    if (fragCount != 0)
    {
        meanFragLen = boost::accumulators::mean(fragLengthAccumulator_);
    }

    return { alleleCount, meanReadLength, meanFragLen, depth };
}

void LocusStatsCalculator::recordReadLen(const GraphAlignment& readAlign)
{
    const graphtools::NodeId firstNode = readAlign.path().getNodeIdByIndex(0);
    if (firstNode == leftFlankId_ || firstNode == rightFlankId_)
    {
        readLengthAccumulator_(readAlign.queryLength());
    }
}

void LocusStatsCalculator::recordFragLen(const GraphAlignment& readAlign, const GraphAlignment& mateAlign)
{
    const auto readStartNode = readAlign.path().getNodeIdByIndex(0);
    const auto mateStartNode = mateAlign.path().getNodeIdByIndex(0);
    const bool matesStartOnLeftFlank = readStartNode == leftFlankId_ && mateStartNode == leftFlankId_;
    const bool matesStartOnRightFlank = readStartNode == rightFlankId_ && mateStartNode == rightFlankId_;

    if (!matesStartOnLeftFlank && !matesStartOnRightFlank)
    {
        return;
    }

    const int readStart = readAlign.path().startPosition();
    const int readEnd = readStart + static_cast<int>(readAlign.queryLength());

    const int mateStart = mateAlign.path().startPosition();
    const int mateEnd = mateStart + static_cast<int>(mateAlign.queryLength());

    if (readEnd < mateEnd)
    {
        fragLengthAccumulator_(mateEnd - readStart);
    }
    else if (mateEnd < readEnd)
    {
        fragLengthAccumulator_(readEnd - mateStart);
    }
}

}
