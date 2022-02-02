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

#include "locus/RepeatAnalyzer.hh"

#include <boost/smart_ptr/make_unique.hpp>

// clang-format off
// Note that spdlog.h must be included before ostr.h
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
// clang-format on

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphutils/SequenceOperations.hh"

#include "alignment/AlignmentFilters.hh"
#include "genotyping/AlignMatrixFiltering.hh"
#include "genotyping/StrGenotyper.hh"

namespace ehunter
{

using boost::make_unique;
using boost::optional;
using graphtools::GraphAlignment;
using graphtools::NodeId;
using graphtools::prettyPrint;
using graphtools::splitStringByDelimiter;
using std::list;
using std::string;
using std::unique_ptr;
using std::vector;

void RepeatAnalyzer::processMates(
    const Read& /*read*/, const GraphAlignment& readAlignment, const Read& /*mate*/,
    const GraphAlignment& mateAlignment)
{
    alignMatrix_.add(readAlignment, mateAlignment);
    alignmentStatsCalculator_.inspect(readAlignment);
    alignmentStatsCalculator_.inspect(mateAlignment);
}

unique_ptr<VariantFindings> RepeatAnalyzer::analyze(const LocusStats& stats)
{
    if (isLowDepth(stats))
    {
        return make_unique<RepeatFindings>(
            CountTable(), CountTable(), CountTable(), stats.alleleCount(), boost::none, GenotypeFilter::kLowDepth);
    }

    auto genotypeFilter = GenotypeFilter();

    const int minBreakpointSpanningReads = stats.alleleCount() == AlleleCount::kTwo
        ? genotyperParams_.minBreakpointSpanningReads
        : (genotyperParams_.minBreakpointSpanningReads / 2);

    GraphVariantAlignmentStats alignmentStats = alignmentStatsCalculator_.getStats();
    if (alignmentStats.numReadsSpanningRightBreakpoint() < minBreakpointSpanningReads
        || alignmentStats.numReadsSpanningLeftBreakpoint() < minBreakpointSpanningReads)
    {
        genotypeFilter = genotypeFilter | GenotypeFilter::kLowDepth;
    }

    if (countOfInrepeatReadPairs_)
    {
        const int maxMotifsInRead = std::ceil(stats.meanReadLength() / static_cast<double>(repeatUnit_.length()));
        strgt::addIrrPairsIfPossibleExpansion(maxMotifsInRead, alignMatrix_, countOfInrepeatReadPairs_);
    }

    auto countsOfSpanningReads = countAligns(StrAlign::Type::kSpanning, alignMatrix_);
    auto countsOfFlankingReads = countAligns(StrAlign::Type::kFlanking, alignMatrix_);
    auto countsOfInrepeatReads = countAligns(StrAlign::Type::kInRepeat, alignMatrix_);

    auto genotype = strgt::genotype(
        stats.alleleCount(), repeatUnit_.length(), stats.meanReadLength(), stats.medianFragLength(), alignMatrix_);

    return make_unique<RepeatFindings>(
        countsOfSpanningReads, countsOfFlankingReads, countsOfInrepeatReads, stats.alleleCount(), genotype,
        genotypeFilter);
}

}
