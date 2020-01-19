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

#include <cassert>
#include <memory>
#include <string>
#include <unordered_map>

#include <boost/optional.hpp>

#include "thirdparty/spdlog/fmt/ostr.h"
#include "thirdparty/spdlog/spdlog.h"

#include "graphalign/GappedAligner.hh"
#include "graphalign/KmerIndex.hh"
#include "graphio/AlignmentWriter.hh"

#include "alignment/SoftclippingAligner.hh"
#include "common/Parameters.hh"
#include "filtering/OrientationPredictor.hh"
#include "reads/Read.hh"
#include "region_analysis/LocusFindings.hh"
#include "region_analysis/VariantAnalyzer.hh"
#include "region_spec/LocusSpecification.hh"
#include "stats/LocusStats.hh"
#include "stats/WeightedPurityCalculator.hh"

namespace ehunter
{

// Regions of the reference genome that can contain reads that originated in a given locus are partitioned into target
// and offtarget regions. Target regions typically consist of the reference region of the locus and possibly other
// highly-similar regions where reads typically misalign. Offtarget regions are regions where certain kinds of
// relevant reads might occasionally misalign and that require special handling (usually for efficiency reasons).
enum class RegionType
{
    kTarget,
    kOfftarget
};

class LocusAnalyzer
{
public:
    LocusAnalyzer(
        const LocusSpecification& locusSpec, HeuristicParameters heuristicParams,
        graphtools::AlignmentWriter& alignmentWriter);

    LocusAnalyzer(const LocusAnalyzer&) = delete;
    LocusAnalyzer& operator=(const LocusAnalyzer&) = delete;
    LocusAnalyzer(LocusAnalyzer&&) = default;

    const std::string& locusId() const { return locusSpec_.locusId(); }
    const LocusSpecification& locusSpec() const { return locusSpec_; }

    void processMates(Read read, boost::optional<Read> mate, RegionType regionType);

    bool checkIfPassesSequenceFilters(const std::string& sequence) const;

    LocusFindings analyze(Sex sampleSex, boost::optional<double> genomeWideDepth);

    bool operator==(const LocusAnalyzer& other) const;

private:
    void processOntargetMates(Read read, boost::optional<Read> mate);
    void processOfftargetMates(Read read, Read mate);
    void runVariantAnalysis(
        const Read& read, const GraphAlignment& readAlignment, const Read& mate, const GraphAlignment& mateAlignment);
    boost::optional<GraphAlignment> alignRead(Read& read) const;

    LocusSpecification locusSpec_;
    HeuristicParameters heuristicParams_;

    graphtools::AlignmentWriter& alignmentWriter_;
    OrientationPredictor orientationPredictor_;
    SoftclippingAligner graphAligner_;

    std::unordered_map<std::string, WeightedPurityCalculator> weightedPurityCalculators;
    LocusStatsCalculator statsCalculator_;

    std::vector<std::unique_ptr<VariantAnalyzer>> variantAnalyzerPtrs_;
    boost::optional<std::string> optionalUnitOfRareRepeat_;

    std::shared_ptr<spdlog::logger> console_;
};

std::vector<std::unique_ptr<LocusAnalyzer>> initializeLocusAnalyzers(
    const RegionCatalog& RegionCatalog, const HeuristicParameters& heuristicParams,
    graphtools::AlignmentWriter& alignmentWriter);

}
