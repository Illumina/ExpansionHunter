//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
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

#include "classification/AlignmentClassifier.hh"

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