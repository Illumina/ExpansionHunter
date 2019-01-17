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

#include "classification/ClassifierOfAlignmentsToVariant.hh"

#include <map>

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"

#include "common/CountTable.hh"
#include "input/GraphBlueprint.hh"
#include "input/RegionGraph.hh"

using graphtools::Graph;
using graphtools::GraphAlignment;
using std::map;
using std::string;

using namespace ehunter;

TEST(InitializingAlignmentClassifier, VariantNodesAreNonconsecutive_ExceptionThrown)
{
    EXPECT_ANY_THROW(ClassifierOfAlignmentsToVariant({}));
    EXPECT_ANY_THROW(ClassifierOfAlignmentsToVariant({ 2, 4 }));
}

TEST(ClassifyingAlignmentsOverIndel, DownstreamAndUpstreamAlignments_Classified)
{
    ClassifierOfAlignmentsToVariant classifier({ 4 });

    //                                          NodeIds =  0  1 2 3  4   5
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("AC(T|G)CT(CA)?TGTGT"));

    GraphAlignment upstreamAlignment = decodeGraphAlignment(1, "0[1M]1[1M]3[2M]", &graph);
    ASSERT_TRUE(checkConsistency(upstreamAlignment, "CTCT"));

    GraphAlignment downstreamAlignment = decodeGraphAlignment(0, "5[4M]", &graph);
    ASSERT_TRUE(checkConsistency(downstreamAlignment, "TGTG"));

    GraphAlignment spanningAlignment = decodeGraphAlignment(1, "0[1M]1[1M]3[2M]4[2M]5[3M]", &graph);
    ASSERT_TRUE(checkConsistency(spanningAlignment, "CTCTCATGT"));

    GraphAlignment bypassingAlignment = decodeGraphAlignment(1, "0[1M]1[1M]3[2M]5[3M]", &graph);
    ASSERT_TRUE(checkConsistency(bypassingAlignment, "CTCTTGT"));

    GraphAlignment upstreamFlankingAlignment = decodeGraphAlignment(1, "0[1M]1[1M]3[2M]4[2M]", &graph);
    ASSERT_TRUE(checkConsistency(upstreamFlankingAlignment, "CTCTCA"));

    GraphAlignment downstreamFlankingAlignment = decodeGraphAlignment(0, "4[2M]5[3M]", &graph);
    ASSERT_TRUE(checkConsistency(downstreamFlankingAlignment, "CATGT"));

    classifier.classify(upstreamAlignment);
    classifier.classify(downstreamAlignment);
    classifier.classify(spanningAlignment);
    classifier.classify(bypassingAlignment);
    classifier.classify(upstreamFlankingAlignment);
    classifier.classify(downstreamFlankingAlignment);

    EXPECT_EQ(CountTable(map<int32_t, int32_t>({ { 4, 1 } })), classifier.countsOfReadsFlankingUpstream());
    EXPECT_EQ(CountTable(map<int32_t, int32_t>({ { 4, 1 } })), classifier.countsOfReadsFlankingDownstream());
    EXPECT_EQ(CountTable(map<int32_t, int32_t>({ { 4, 1 } })), classifier.countsOfSpanningReads());
    EXPECT_EQ(1, classifier.numBypassingReads());
}
