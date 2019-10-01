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

#include "classification/StrAlignmentClassifier.hh"

#include <list>
#include <memory>
#include <string>

#include "gtest/gtest.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"

using graphtools::Graph;
using graphtools::makeStrGraph;
using graphtools::NodeId;
using std::list;
using std::string;

using namespace ehunter;

TEST(InitializingAlignmentClassifier, SingleUnitRepeatGraph_RepeatFlanksDetermined)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    StrAlignmentClassifier classifier(graph, 1);

    const std::set<NodeId> expectedLeftFlankIds = { 0 };
    const std::set<NodeId> expectedRightFlankIds = { 2 };

    EXPECT_EQ(expectedLeftFlankIds, classifier.leftFlankNodeIds());
    EXPECT_EQ(expectedRightFlankIds, classifier.rightFlankNodeIds());
}

TEST(AlignmentClassificaton, SpanningAlignment_Classified)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    StrAlignmentClassifier classifier(graph, 1);

    { //                     FFRRRRRRFF
        const string read = "CCCCGCCGAT";
        GraphAlignment alignment = decodeGraphAlignment(4, "0[2M]1[3M]1[3M]2[2M]", &graph);

        int expectedScore = scoreAlignment(alignment);
        StrAlignment expectedSummary(2, StrAlignment::Type::kSpanning, expectedScore, 0);
        EXPECT_EQ(expectedSummary, *classifier.classify(alignment));
    }

    { //                     FFFF
        const string read = "CCAT";
        GraphAlignment alignment = decodeGraphAlignment(4, "0[2M]2[2M]", &graph);

        int expectedScore = scoreAlignment(alignment);
        StrAlignment expectedSummary(0, StrAlignment::Type::kSpanning, expectedScore, 0);
        EXPECT_EQ(expectedSummary, *classifier.classify(alignment));
    }
}

TEST(AlignmentClassificaton, FlankingAlignment_Classified)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    StrAlignmentClassifier classifier(graph, 1);

    { //                  FFFFRRR
        const string read = "AACCCCG";
        GraphAlignment alignment = decodeGraphAlignment(2, "0[4M]1[3M]", &graph);

        int expectedScore = scoreAlignment(alignment);
        StrAlignment expectedSummary(1, StrAlignment::Type::kFlanking, expectedScore, 0);
        EXPECT_EQ(expectedSummary, *classifier.classify(alignment));
    }

    { //                  RRRFFF
        const string read = "CCGATT";
        GraphAlignment alignment = decodeGraphAlignment(0, "1[3M]2[3M]", &graph);

        int expectedScore = scoreAlignment(alignment);
        StrAlignment expectedSummary(1, StrAlignment::Type::kFlanking, expectedScore, 0);
        EXPECT_EQ(expectedSummary, *classifier.classify(alignment));
    }
}

TEST(AlignmentClassificaton, RepeatAlignment_Classified)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    StrAlignmentClassifier classifier(graph, 1);

    { //                  RRRRRRRR
        const string read = "CCGCCGCC";
        GraphAlignment alignment = decodeGraphAlignment(0, "1[3M]1[3M]1[2M]", &graph);

        int expectedScore = scoreAlignment(alignment);
        StrAlignment expectedSummary(2, StrAlignment::Type::kInrepeat, expectedScore, 0);
        EXPECT_EQ(expectedSummary, *classifier.classify(alignment));
    }

    { //                  RRRRRRRR
        const string read = "CGCCGCCG";
        GraphAlignment alignment = decodeGraphAlignment(1, "1[2M]1[3M]1[3M]", &graph);

        int expectedScore = scoreAlignment(alignment);
        StrAlignment expectedSummary(2, StrAlignment::Type::kInrepeat, expectedScore, 0);
        EXPECT_EQ(expectedSummary, classifier.classify(alignment));
    }
}

TEST(AlignmentClassificaton, OutsideRepeatAlignment_Classified)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    StrAlignmentClassifier classifier(graph, 1);

    { //                     FFFFF
        const string read = "AAAAC";
        GraphAlignment alignment = decodeGraphAlignment(0, "0[5M]", &graph);
        EXPECT_FALSE(classifier.classify(alignment));
    }

    { //                     FFF
        const string read = "TTT";
        GraphAlignment alignment = decodeGraphAlignment(1, "2[3M]", &graph);
        EXPECT_FALSE(classifier.classify(alignment));
    }
}
