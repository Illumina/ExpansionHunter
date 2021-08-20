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

#include "gtest/gtest.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"

#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using graphtools::Graph;
using graphtools::GraphAlignment;

using namespace ehunter;

TEST(DISABLE_RegionAnalysis, ShortSingleUnitRepeat_Genotyped)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    GraphAlignment alignmentA1 = decodeGraphAlignment(3, "0[3M]1[1M]1[1M]1[1M]2[4M]", &graph);
    GraphAlignment alignmentA2 = decodeGraphAlignment(4, "0[2M]1[1M]1[1M]1[1M]2[5M]", &graph);
    GraphAlignment alignmentB1 = decodeGraphAlignment(3, "0[3M]1[1M]2[4M]", &graph);
    GraphAlignment alignmentB2 = decodeGraphAlignment(4, "0[2M]1[1M]2[5M]", &graph);
    GraphAlignment alignmentC1 = decodeGraphAlignment(2, "0[4M]1[1M]", &graph);
    GraphAlignment alignmentC2 = decodeGraphAlignment(0, "1[1M]1[1M]2[5M]", &graph);

    // AlleleCount expectedAlleleCount = AlleleCount::kTwo;
    // graphtools::NodeId repeatNodeId = 1;
    // const int32_t maxNumUnitsInRead = 10;
    // const double haplotypeDepth = 5.0;

    // RepeatAnalyzer repeatAnalyzer("repeat0", expectedAlleleCount, graph, repeatNodeId);
    // repeatAnalyzer.processMates({ alignmentA1 }, { alignmentA2 });
    // repeatAnalyzer.processMates({ alignmentB1 }, { alignmentB2 });
    // repeatAnalyzer.processMates({ alignmentC1 }, { alignmentC2 });
    // RepeatFindings repeatFindings = repeatAnalyzer.analyzeCommonRepeat();

    // RepeatGenotype repeatGenotype(1, { 1, 3 });
    // RepeatFindings expectedFindings(
    //    CountTable({ { 1, 1 }, { 2, 1 } }), CountTable({ { 1, 2 }, { 3, 2 } }), repeatGenotype);
    // ASSERT_EQ(expectedFindings, repeatFindings);
}
