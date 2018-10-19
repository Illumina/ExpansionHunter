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

#include "region_analysis/RepeatAnalyzer.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"

#include "region_spec/region_graph.h"

using graphtools::Graph;
using graphtools::GraphAlignment;

TEST(DISABLE_RegionAnalysis, ShortSingleUnitRepeat_Genotyped)
{
    Graph graph = makeRegionGraph("ATTCGA", "C", "ATGTCG");

    GraphAlignment alignmentA1 = decodeGraphAlignment(3, "0[3M]1[1M]1[1M]1[1M]2[4M]", &graph);
    GraphAlignment alignmentA2 = decodeGraphAlignment(4, "0[2M]1[1M]1[1M]1[1M]2[5M]", &graph);
    GraphAlignment alignmentB1 = decodeGraphAlignment(3, "0[3M]1[1M]2[4M]", &graph);
    GraphAlignment alignmentB2 = decodeGraphAlignment(4, "0[2M]1[1M]2[5M]", &graph);
    GraphAlignment alignmentC1 = decodeGraphAlignment(2, "0[4M]1[1M]", &graph);
    GraphAlignment alignmentC2 = decodeGraphAlignment(0, "1[1M]1[1M]2[5M]", &graph);

    AlleleCount expectedAlleleCount = AlleleCount::kTwo;
    graphtools::NodeId repeatNodeId = 1;
    const int32_t maxNumUnitsInRead = 10;
    const double haplotypeDepth = 5.0;

    RepeatAnalyzer repeatAnalyzer(
        "repeat0", graph, expectedAlleleCount, repeatNodeId, haplotypeDepth, maxNumUnitsInRead);
    // repeatAnalyzer.processMates({ alignmentA1 }, { alignmentA2 });
    // repeatAnalyzer.processMates({ alignmentB1 }, { alignmentB2 });
    // repeatAnalyzer.processMates({ alignmentC1 }, { alignmentC2 });
    // RepeatFindings repeatFindings = repeatAnalyzer.analyzeCommonRepeat();

    // RepeatGenotype repeatGenotype(1, { 1, 3 });
    // RepeatFindings expectedFindings(
    //    CountTable({ { 1, 1 }, { 2, 1 } }), CountTable({ { 1, 2 }, { 3, 2 } }), repeatGenotype);
    // ASSERT_EQ(expectedFindings, repeatFindings);
}
