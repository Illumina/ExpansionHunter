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

#include "classification/SmallVariantAlignmentClassifier.hh"

#include <map>

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"

#include "common/CountTable.hh"
#include "locus_spec/GraphBlueprint.hh"
#include "locus_spec/RegionGraph.hh"

using boost::optional;
using graphtools::Graph;
using graphtools::GraphAlignment;
using std::map;
using std::string;

using namespace ehunter;

TEST(InitializingAlignmentClassifier, VariantNodesAreNonconsecutive_ExceptionThrown)
{
    EXPECT_ANY_THROW(SmallVariantAlignmentClassifier({}));
    EXPECT_ANY_THROW(SmallVariantAlignmentClassifier({ 2, 4 }));
}

TEST(ClassifyingAlignmentsOverIndel, DownstreamAndUpstreamAlignments_Classified)
{
    //                                          NodeIds =  0  1 2 3  4   5
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("AC(T|G)CT(CA)?TGTGT"));
    SmallVariantAlignmentClassifier classifier({ 4 });

    {
        GraphAlignment alignment = decodeGraphAlignment(1, "0[1M]1[1M]3[2M]", &graph);
        string read = "CTCT";
        ASSERT_TRUE(checkConsistency(alignment, read));

        ReadSummaryForSmallVariant expectedRead(read.length());
        EXPECT_EQ(expectedRead, classifier.classifyRead(read, { alignment }));
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(0, "5[4M]", &graph);
        string read = "TGTG";
        ASSERT_TRUE(checkConsistency(alignment, read));

        ReadSummaryForSmallVariant expectedRead(read.length());
        EXPECT_EQ(expectedRead, classifier.classifyRead(read, { alignment }));
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(1, "0[1M]1[1M]3[2M]4[2M]5[3M]", &graph);
        string read = "CTCTCATGT";
        ASSERT_TRUE(checkConsistency(alignment, read));

        ReadSummaryForSmallVariant expectedSummary(read.length());
        int expectedScore = scoreAlignment(alignment);
        expectedSummary.addAlignment(SmallVariantAlignment(4, SmallVariantAlignment::Type::kSpanning, expectedScore));
        EXPECT_EQ(expectedSummary, classifier.classifyRead(read, { alignment }));
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(1, "0[1M]1[1M]3[2M]5[3M]", &graph);
        string read = "CTCTTGT";
        ASSERT_TRUE(checkConsistency(alignment, read));

        ReadSummaryForSmallVariant expectedSummary(read.length());
        int expectedNode = SmallVariantAlignmentClassifier::kInvalidNodeId;
        int expectedScore = scoreAlignment(alignment);
        expectedSummary.addAlignment(
            SmallVariantAlignment(expectedNode, SmallVariantAlignment::Type::kSpanning, expectedScore));
        EXPECT_EQ(expectedSummary, classifier.classifyRead(read, { alignment }));
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(1, "0[1M]1[1M]3[2M]4[2M]", &graph);
        string read = "CTCTCA";
        ASSERT_TRUE(checkConsistency(alignment, read));

        ReadSummaryForSmallVariant expectedSummary(read.length());
        int expectedScore = scoreAlignment(alignment);
        expectedSummary.addAlignment(
            SmallVariantAlignment(4, SmallVariantAlignment::Type::kUpstreamFlanking, expectedScore));
        EXPECT_EQ(expectedSummary, classifier.classifyRead(read, { alignment }));
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(0, "4[2M]5[3M]", &graph);
        string read = "CATGT";
        ASSERT_TRUE(checkConsistency(alignment, read));

        ReadSummaryForSmallVariant expectedSummary(read.length());
        int expectedScore = scoreAlignment(alignment);
        expectedSummary.addAlignment(
            SmallVariantAlignment(4, SmallVariantAlignment::Type::kDownstreamFlanking, expectedScore));
        EXPECT_EQ(expectedSummary, classifier.classifyRead(read, { alignment }));
    }
}
