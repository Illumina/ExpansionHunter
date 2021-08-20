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

#include "alignment/AlignmentClassifier.hh"

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
    RepeatAlignmentClassifier alignment_classifier(graph, 1);

    const std::set<NodeId> expected_left_flank_ids = { 0 };
    const std::set<NodeId> expected_right_flank_ids = { 2 };

    EXPECT_EQ(expected_left_flank_ids, alignment_classifier.leftFlankNodeIds());
    EXPECT_EQ(expected_right_flank_ids, alignment_classifier.rightFlankNodeIds());
}

TEST(AlignmentClassificaton, SpanningAlignment_Classified)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    RepeatAlignmentClassifier alignment_classifier(graph, 1);

    { //                  FFRRRRRRFF
        const string read = "CCCCGCCGAT";
        GraphAlignment alignment = decodeGraphAlignment(4, "0[2M]1[3M]1[3M]2[2M]", &graph);

        EXPECT_EQ(AlignmentType::kSpansRepeat, alignment_classifier.Classify(alignment));
    }

    { //                  FFFF
        const string read = "CCAT";
        GraphAlignment alignment = decodeGraphAlignment(4, "0[2M]2[2M]", &graph);

        RepeatAlignmentClassifier alignment_classifier(graph, 1);
        EXPECT_EQ(AlignmentType::kSpansRepeat, alignment_classifier.Classify(alignment));
    }
}

TEST(AlignmentClassificaton, FlankingAlignment_Classified)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    RepeatAlignmentClassifier alignment_classifier(graph, 1);

    { //                  FFFFRRR
        const string read = "AACCCCG";
        GraphAlignment alignment = decodeGraphAlignment(2, "0[4M]1[3M]", &graph);

        EXPECT_EQ(AlignmentType::kFlanksRepeat, alignment_classifier.Classify(alignment));
    }

    { //                  RRRFFF
        const string read = "CCGATT";
        GraphAlignment alignment = decodeGraphAlignment(0, "1[3M]2[3M]", &graph);

        RepeatAlignmentClassifier alignment_classifier(graph, 1);
        EXPECT_EQ(AlignmentType::kFlanksRepeat, alignment_classifier.Classify(alignment));
    }
}

TEST(AlignmentClassificaton, RepeatAlignment_Classified)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    RepeatAlignmentClassifier alignment_classifier(graph, 1);

    { //                  RRRRRRRR
        const string read = "CCGCCGCC";
        GraphAlignment alignment = decodeGraphAlignment(0, "1[3M]1[3M]1[2M]", &graph);

        EXPECT_EQ(AlignmentType::kInsideRepeat, alignment_classifier.Classify(alignment));
    }

    { //                  RRRRRRRR
        const string read = "CGCCGCCG";
        GraphAlignment alignment = decodeGraphAlignment(1, "1[2M]1[3M]1[3M]", &graph);

        RepeatAlignmentClassifier alignment_classifier(graph, 1);
        EXPECT_EQ(AlignmentType::kInsideRepeat, alignment_classifier.Classify(alignment));
    }
}

TEST(AlignmentClassificaton, OutsideRepeatAlignment_Classified)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    RepeatAlignmentClassifier alignment_classifier(graph, 1);

    { //                  FFFFF
        const string read = "AAAAC";
        GraphAlignment alignment = decodeGraphAlignment(0, "0[5M]", &graph);

        EXPECT_EQ(AlignmentType::kOutsideRepeat, alignment_classifier.Classify(alignment));
    }

    { //                  FFF
        const string read = "TTT";
        GraphAlignment alignment = decodeGraphAlignment(1, "2[3M]", &graph);

        RepeatAlignmentClassifier alignment_classifier(graph, 1);
        EXPECT_EQ(AlignmentType::kOutsideRepeat, alignment_classifier.Classify(alignment));
    }
}

TEST(ObtainingCanonicalAlignment, FlankingAndRepeatRead_ClassifiedAsRepeat)
{
    Graph graph = makeStrGraph("AAAACG", "CCG", "ATTT");
    RepeatAlignmentClassifier alignment_classifier(graph, 1);

    //                   FFFFFFFF
    const string read = "CGCCGCCG";
    const GraphAlignment flanking_alignment = decodeGraphAlignment(4, "0[2M]1[3M]1[3M]", &graph);
    const GraphAlignment irr_alignment = decodeGraphAlignment(1, "1[2M]1[3M]1[3M]", &graph);

    const list<GraphAlignment> alignments = { irr_alignment, flanking_alignment };

    EXPECT_EQ(irr_alignment, alignment_classifier.GetCanonicalAlignment(alignments));
}

TEST(ObtainingCanonicalAlignment, FlankingAndSpanningRead_ClassifiedAsFlanking)
{
    Graph graph = makeStrGraph("AAAACG", "CCG", "ATTT");
    RepeatAlignmentClassifier alignment_classifier(graph, 1);

    //                   FFFFFFFFFF
    const string read = "CGCCGCCGAT";
    const GraphAlignment spanning_Alignment = decodeGraphAlignment(4, "0[2M]1[3M]1[3M]2[2M]", &graph);
    const GraphAlignment flanking_Alignment = decodeGraphAlignment(1, "1[2M]1[3M]1[3M]2[2M]", &graph);

    const list<GraphAlignment> Alignments = { spanning_Alignment, flanking_Alignment };

    EXPECT_EQ(flanking_Alignment, alignment_classifier.GetCanonicalAlignment(Alignments));
}
