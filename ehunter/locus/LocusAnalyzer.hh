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

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "graphalign/GappedAligner.hh"
#include "graphio/AlignmentWriter.hh"

#include "core/GenomicRegion.hh"
#include "core/LocusStats.hh"
#include "core/Read.hh"
#include "locus/AlignmentBuffer.hh"
#include "locus/IrrPairFinder.hh"
#include "locus/LocusAligner.hh"
#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"
#include "locus/VariantAnalyzer.hh"

namespace ehunter
{
namespace locus
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

using AlignWriterPtr = std::shared_ptr<graphtools::AlignmentWriter>;

class LocusAnalyzer
{
public:
    using Align = graphtools::GraphAlignment;
    using Node = graphtools::NodeId;

    LocusAnalyzer(LocusSpecification locusSpec, const HeuristicParameters& params, AlignWriterPtr writer);

    const std::string& locusId() const { return locusSpec_.locusId(); }
    const LocusSpecification& locusSpec() const { return locusSpec_; }

    void processMates(Read& read, Read* mate, RegionType regionType, graphtools::AlignerSelector& alignerSelector);
    LocusFindings analyze(Sex sampleSex, boost::optional<double> genomeWideDepth);

    const boost::optional<IrrPairFinder>& irrPairFinder() const { return irrPairFinder_; }
    void addIrrPairFinder(std::string motif);

    void addRepeatAnalyzer(std::string variantId, Node nodeId);
    void addSmallVariantAnalyzer(
        std::string variantId, VariantSubtype subtype, std::vector<Node> nodes, boost::optional<Node> refNode);

private:
    void processOntargetMates(Read& read, Read* mate, graphtools::AlignerSelector& alignerSelector);
    void processOfftargetMates(const Read& read, const Read& mate);
    void runVariantAnalysis(const Read& read, const Align& readAlign, const Read& mate, const Align& mateAlign);

    LocusSpecification locusSpec_;

    // Read alignments are optionally buffered for custom additional analysis at certain loci
    std::shared_ptr<locus::AlignmentBuffer> alignmentBuffer_;

    LocusAligner aligner_;
    LocusStatsCalculator statsCalc_;
    boost::optional<IrrPairFinder> irrPairFinder_;
    std::vector<std::unique_ptr<VariantAnalyzer>> variantAnalyzers_;
};

}
}
