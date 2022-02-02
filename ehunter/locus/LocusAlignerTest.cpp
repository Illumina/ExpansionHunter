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

#include "locus/LocusAligner.hh"

#include "gmock/gmock.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphio/AlignmentWriter.hh"

#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using namespace ehunter;
using namespace locus;

using boost::optional;
using graphtools::BlankAlignmentWriter;
using graphtools::decodeGraphAlignment;
using graphtools::GraphAlignment;

LocusAligner makeStrAligner(graphtools::Graph* graph)
{
    HeuristicParameters params(1000, 10, 20, true, graphtools::AlignerType::DAG_ALIGNER, 4, 0, 0, 4, 1);
    auto writer = std::make_shared<BlankAlignmentWriter>();
    return { "str", graph, params, writer, {} };
}

TEST(AligningReads, ReadPairFromSameStrand_Aligned)
{
    Read read(ReadId("frag1", MateNumber::kFirstMate), "ATTACC", true);
    Read mate(ReadId("frag1", MateNumber::kSecondMate), "GGCGGC", true);

    auto graph = makeRegionGraph(decodeFeaturesFromRegex("ATATTA(C)*GGCGGC"));
    auto aligner = makeStrAligner(&graph);
    graphtools::AlignerSelector selector(graphtools::AlignerType::DAG_ALIGNER);
    auto alignedPair = aligner.align(read, &mate, selector);

    GraphAlignment expectedReadAlign = decodeGraphAlignment(2, "0[4M]1[1M]1[1M]", &graph);
    GraphAlignment expectedMateAlign = decodeGraphAlignment(0, "2[6M]", &graph);

    ASSERT_TRUE(alignedPair.first && alignedPair.second); // Both reads must be aligned
    ASSERT_EQ(expectedReadAlign, *alignedPair.first);
    ASSERT_EQ(expectedMateAlign, *alignedPair.second);
}

TEST(AligningReads, ReadPairFromOppositeStrand_Aligned)
{
    Read read(ReadId("frag1", MateNumber::kFirstMate), "GGTAAT", true);
    Read mate(ReadId("frag1", MateNumber::kSecondMate), "GCCGCC", false);

    auto graph = makeRegionGraph(decodeFeaturesFromRegex("ATATTA(C)*GGCGGC"));
    auto aligner = makeStrAligner(&graph);
    graphtools::AlignerSelector selector(graphtools::AlignerType::DAG_ALIGNER);
    auto alignedPair = aligner.align(read, &mate, selector);

    GraphAlignment expectedReadAlign = decodeGraphAlignment(2, "0[4M]1[1M]1[1M]", &graph);
    GraphAlignment expectedMateAlign = decodeGraphAlignment(0, "2[6M]", &graph);

    ASSERT_TRUE(alignedPair.first && alignedPair.second); // Both reads must be aligned
    ASSERT_EQ(expectedReadAlign, *alignedPair.first);
    ASSERT_EQ(expectedMateAlign, *alignedPair.second);
    ASSERT_FALSE(read.isReversed());
    ASSERT_TRUE(mate.isReversed());
}

TEST(AligningReads, NotLocallyPlasedReads_NotAligned)
{
    Read read(ReadId("frag1", MateNumber::kFirstMate), "TACCC", true);
    Read mate(ReadId("frag1", MateNumber::kSecondMate), "CCCGG", false);

    auto graph = makeRegionGraph(decodeFeaturesFromRegex("ATATTA(C)*GGCGGC"));
    auto aligner = makeStrAligner(&graph);
    graphtools::AlignerSelector selector(graphtools::AlignerType::DAG_ALIGNER);
    auto alignedPair = aligner.align(read, &mate, selector);

    ASSERT_FALSE(alignedPair.first);
    ASSERT_FALSE(alignedPair.second);
}
