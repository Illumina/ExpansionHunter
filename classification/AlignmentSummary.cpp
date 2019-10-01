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

#include "classification/AlignmentSummary.hh"

#include <stdexcept>

#include "graphalign/LinearAlignmentOperations.hh"

namespace ehunter
{

using graphtools::GraphAlignment;
using std::ostream;

ostream& operator<<(ostream& os, const StrAlignment::Type& alignmentType)
{
    switch (alignmentType)
    {
    case StrAlignment::Type::kSpanning:
        os << "kSpanning";
        break;
    case StrAlignment::Type::kFlanking:
        os << "kFlanking";
        break;
    case StrAlignment::Type::kInrepeat:
        os << "kInrepeat";
        break;
    default:
        throw std::logic_error("Encountered unknown StrAlignment::Type");
    }
    return os;
}

ostream& operator<<(ostream& os, const StrAlignment& alignmentSummary)
{
    os << "StrAlignment(" << alignmentSummary.numUnits() << ", " << alignmentSummary.type() << ", "
       << alignmentSummary.score() << ")";
    return os;
}

ostream& operator<<(ostream& os, const SmallVariantAlignment::Type& alignmentType)
{
    switch (alignmentType)
    {
    case SmallVariantAlignment::Type::kUpstreamFlanking:
        os << "kUpstreamFlanking";
        break;
    case SmallVariantAlignment::Type::kSpanning:
        os << "kSpanning";
        break;
    case SmallVariantAlignment::Type::kDownstreamFlanking:
        os << "kDownstreamFlanking";
        break;
    default:
        throw std::logic_error("Encountered unknown SmallVariantAlignment::Type");
    }
    return os;
}

ostream& operator<<(ostream& os, const SmallVariantAlignment& alignmentSummary)
{
    os << "SmallVariantAlignment(" << alignmentSummary.nodeId() << ", " << alignmentSummary.type() << ", "
       << alignmentSummary.score() << ")";
    return os;
}

int scoreAlignment(const GraphAlignment& alignment, LinearAlignmentParameters parameters)
{
    int score = 0;

    for (int nodeIndex = 0; nodeIndex != static_cast<int>(alignment.size()); ++nodeIndex)
    {
        score += graphtools::scoreAlignment(
            alignment[nodeIndex], parameters.matchScore, parameters.mismatchScore, parameters.gapOpenScore);
    }

    return score;
}

}
