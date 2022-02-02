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

#include "locus/LocusAnalyzer.hh"

#include "gtest/gtest.h"

#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using namespace ehunter;
using namespace locus;

using graphtools::AlignerType;
using std::vector;

LocusSpecification buildStrSpec(const std::string& structure)
{
    auto graph = makeRegionGraph(decodeFeaturesFromRegex(structure));
    vector<GenomicRegion> referenceRegions = { GenomicRegion(1, 1, 2) };

    NodeToRegionAssociation dummyAssociation;
    GenotyperParameters params(10);
    const bool useRFC1MotifAnalysis(false);
    LocusSpecification locusSpec(
        "region", ChromType::kAutosome, referenceRegions, graph, dummyAssociation, params, useRFC1MotifAnalysis);
    VariantClassification classification(VariantType::kRepeat, VariantSubtype::kCommonRepeat);
    locusSpec.addVariantSpecification("repeat", classification, GenomicRegion(1, 1, 2), { 1 }, 1);
    return locusSpec;
}

TEST(CreatingLocusAnalyzer, TypicalParameters_Created)
{
    auto locusSpec = buildStrSpec("ATTCGA(C)*ATGTCG");

    HeuristicParameters heuristicParams(1000, 10, 20, true, AlignerType::DAG_ALIGNER, 4, 1, 5, 4, 1);
    auto writer = std::make_shared<graphtools::BlankAlignmentWriter>();

    graphtools::AlignerSelector selector(heuristicParams.alignerType());
    LocusAnalyzer locusAnalyzer(locusSpec, heuristicParams, writer);

    Read read1(ReadId("read1", MateNumber::kFirstMate), "CGACCCATGT", true);
    Read mate1(ReadId("read1", MateNumber::kSecondMate), "GACCCATGTC", true);
    locusAnalyzer.processMates(read1, &mate1, RegionType::kTarget, selector);

    Read read2(ReadId("read2", MateNumber::kFirstMate), "CGACATGT", true);
    Read mate2(ReadId("read2", MateNumber::kSecondMate), "GACATGTC", true);
    locusAnalyzer.processMates(read2, &mate2, RegionType::kTarget, selector);

    LocusFindings locusFindings = locusAnalyzer.analyze(Sex::kFemale, boost::none);
    RepeatFindings observed = *dynamic_cast<RepeatFindings*>(locusFindings.findingsForEachVariant["repeat"].get());

    CountTable spanningCounts({ { 1, 2 }, { 3, 2 } });
    RepeatGenotype genotype(1, { 1, 3 });
    RepeatFindings repeatFindings(spanningCounts, {}, {}, AlleleCount::kTwo, genotype, {});

    ASSERT_EQ(repeatFindings, observed);
}
