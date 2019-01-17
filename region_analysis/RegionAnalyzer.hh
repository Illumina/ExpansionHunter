//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
#include "region_analysis/VariantAnalyzer.hh"
#include "region_analysis/VariantFindings.hh"
#include "region_spec/LocusSpecification.hh"
#include "stats/WeightedPurityCalculator.hh"

namespace ehunter
{

class RegionAnalyzer
{
public:
    RegionAnalyzer(
        const LocusSpecification& regionSpec, HeuristicParameters heuristicParams,
        graphtools::AlignmentWriter& alignmentWriter);

    RegionAnalyzer(const RegionAnalyzer&) = delete;
    RegionAnalyzer& operator=(const RegionAnalyzer&) = delete;
    RegionAnalyzer(RegionAnalyzer&&) = default;
    RegionAnalyzer& operator=(RegionAnalyzer&&) = default;

    const std::string& regionId() const { return regionSpec_.regionId(); }
    const LocusSpecification& regionSpec() const { return regionSpec_; }

    void processMates(Read read, Read mate);
    void processOfftargetMates(Read read1, Read read2);
    bool checkIfPassesSequenceFilters(const std::string& sequence) const;

    RegionFindings analyze(const SampleParameters& params);

    bool operator==(const RegionAnalyzer& other) const;

private:
    boost::optional<GraphAlignment> alignRead(Read& read) const;

    LocusSpecification regionSpec_;
    HeuristicParameters heuristicParams_;

    graphtools::AlignmentWriter& alignmentWriter_;
    OrientationPredictor orientationPredictor_;
    SoftclippingAligner graphAligner_;

    std::unordered_map<std::string, WeightedPurityCalculator> weightedPurityCalculators;

    std::vector<std::unique_ptr<VariantAnalyzer>> variantAnalyzerPtrs_;
    boost::optional<std::string> optionalUnitOfRareRepeat_;

    std::shared_ptr<spdlog::logger> console_;
};

std::vector<std::unique_ptr<RegionAnalyzer>> initializeRegionAnalyzers(
    const RegionCatalog& RegionCatalog, const HeuristicParameters& heuristicParams,
    graphtools::AlignmentWriter& alignmentWriter);

}
