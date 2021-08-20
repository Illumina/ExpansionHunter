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

#include "alignment//ClassifierOfAlignmentsToVariant.hh"

#include <map>

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"

#include "core/CountTable.hh"
#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

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
