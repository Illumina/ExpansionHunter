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

#include "graph_components/ReadClassifier.hh"

using std::vector;

namespace ehunter
{

ReadClassifier::ReadClassifier(vector<GenomicRegion> targetRegions)
    : targetRegions_(std::move(targetRegions))
{
}

RegionProximity ReadClassifier::classify(const MappedRead& read) const
{
    RegionProximity proximity = RegionProximity::kFar;
    const int64_t readEnd = read.approximateEnd();
    for (const auto& region : targetRegions_)
    {
        if (read.contigIndex() != region.contigIndex())
        {
            continue;
        }

        if (region.start() <= read.pos() && readEnd <= region.end())
        {
            return RegionProximity::kInside;
        }

        if (region.start() - kMinOfftargetDistance_ <= read.pos() && readEnd <= region.end() + kMinOfftargetDistance_)
        {
            proximity = RegionProximity::kOverlapsOrNear;
        }
    }

    return proximity;
}

RegionProximity ReadClassifier::classify(const MappedRead& read, const MappedRead& mate) const
{
    RegionProximity readProximity = classify(read);
    RegionProximity mateProximity = classify(mate);

    if (readProximity == RegionProximity::kInside || mateProximity == RegionProximity::kInside)
    {
        return RegionProximity::kInside;
    }

    if (readProximity == RegionProximity::kOverlapsOrNear || mateProximity == RegionProximity::kOverlapsOrNear)
    {
        return RegionProximity::kOverlapsOrNear;
    }

    return RegionProximity::kFar;
}

std::ostream& operator<<(std::ostream& out, RegionProximity type)
{
    switch (type)
    {
    case RegionProximity::kInside:
        out << "RegionProximity::kInside";
        break;
    case RegionProximity::kOverlapsOrNear:
        out << "RegionProximity::kOverlapsOrNear";
        break;
    case RegionProximity::kFar:
        out << "RegionProximity::kFar";
        break;
    }

    return out;
}

}