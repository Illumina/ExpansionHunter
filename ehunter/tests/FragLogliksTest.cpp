//
// ExpansionHunter
// Copyright 2016-2020 Illumina, Inc.
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

#include "genotyping/FragLogliks.hh"

#include "gmock/gmock.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"

#include "genotyping/AlignMatrix.hh"
#include "genotyping/RepeatGenotype.hh"
#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using namespace ehunter;
using namespace strgt;
using graphtools::Graph;
using graphtools::GraphAlignment;
using std::vector;
using testing::DoubleNear;
using testing::Pointwise;

TEST(ReadLogliks, TypicalReads_Computed)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    AlignMatrix alignMatrix(1);
    GraphAlignment mate = decodeGraphAlignment(0, "0[6M]", &graph);
    alignMatrix.add(decodeGraphAlignment(3, "0[3M]1[1M]1[1M]2[4M]", &graph), mate);
    alignMatrix.add(decodeGraphAlignment(3, "0[3M]1[1M]1[1M]", &graph), mate);
    alignMatrix.add(decodeGraphAlignment(0, "1[1M]1[1M]1[1M]2[4M]", &graph), mate);
    alignMatrix.add(mate, decodeGraphAlignment(0, "1[1M]1[1M]1[1M]1[1M]", &graph));
    ASSERT_EQ(alignMatrix.numReads(), 8);

    int motifLen = 1;
    int readLen = 8;
    int fragLen = 20;

    FragLogliks logliks(motifLen, readLen, fragLen, &alignMatrix);

    EXPECT_THAT(logliks.getLoglik(0, 0), DoubleNear(-14.45, 0.1));
    EXPECT_THAT(logliks.getLoglik(0, 1), DoubleNear(-13.23, 0.1));
    EXPECT_THAT(logliks.getLoglik(0, 2), DoubleNear(-8.07, 0.1));
    EXPECT_THAT(logliks.getLoglik(0, 3), DoubleNear(-12.10, 0.1));
    EXPECT_THAT(logliks.getLoglik(0, 3), DoubleNear(-12.10, 0.1)); // Test return of cached values

    EXPECT_THAT(logliks.getLoglik(1, 0), DoubleNear(-15.76, 0.1));
    EXPECT_THAT(logliks.getLoglik(1, 1), DoubleNear(-14.55, 0.1));
    EXPECT_THAT(logliks.getLoglik(1, 2), DoubleNear(-13.32, 0.1));
    EXPECT_THAT(logliks.getLoglik(1, 3), DoubleNear(-13.41, 0.1));

    EXPECT_THAT(logliks.getLoglik(2, 0), DoubleNear(-14.45, 0.1));
    EXPECT_THAT(logliks.getLoglik(2, 1), DoubleNear(-13.23, 0.1));
    EXPECT_THAT(logliks.getLoglik(2, 2), DoubleNear(-12.01, 0.1));
    EXPECT_THAT(logliks.getLoglik(2, 3), DoubleNear(-10.79, 0.1));

    EXPECT_THAT(logliks.getLoglik(3, 0), DoubleNear(-19.70, 0.1));
    EXPECT_THAT(logliks.getLoglik(3, 1), DoubleNear(-18.48, 0.1));
    EXPECT_THAT(logliks.getLoglik(3, 2), DoubleNear(-17.26, 0.1));
    EXPECT_THAT(logliks.getLoglik(3, 3), DoubleNear(-16.03, 0.1));
}
