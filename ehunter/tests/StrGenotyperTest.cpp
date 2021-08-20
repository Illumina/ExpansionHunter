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

#include "genotyping/StrGenotyper.hh"

#include <unordered_set>

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
using std::unordered_set;
using std::vector;

TEST(StrAlleleCandidates, TypicalAlignments_Computed)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(CAG)*ATGTCG"));

    AlignMatrix alignMatrix(1);
    int readLen = 24;
    int motifLen = 3;
    GraphAlignment mate = decodeGraphAlignment(0, "0[6M]", &graph);
    alignMatrix.add(decodeGraphAlignment(3, "0[3M]1[3M]1[3M]2[4M]", &graph), mate);
    alignMatrix.add(decodeGraphAlignment(3, "0[3M]1[3M]1[3M]1[3M]2[2M]", &graph), mate);
    ASSERT_EQ(getAlleleCandidates(readLen, motifLen, alignMatrix), unordered_set<int>({ 2, 3 }));

    alignMatrix.add(decodeGraphAlignment(0, "1[3M]1[3M]1[3M]1[3M]1[3M]2[2M]", &graph), mate);
    alignMatrix.add(decodeGraphAlignment(3, "0[3M]1[3M]1[3M]1[3M]", &graph), mate);
    ASSERT_EQ(getAlleleCandidates(readLen, motifLen, alignMatrix), unordered_set<int>({ 2, 3, 5 }));

    /*
    alignMatrix.add(decodeGraphAlignment(0, "1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]", &graph), mate);
    alignMatrix.add(decodeGraphAlignment(0, "1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]", &graph), mate);
    alignMatrix.add(decodeGraphAlignment(0, "1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]", &graph), mate);
    ASSERT_EQ(getStrAlleleCandidates(readLen, motifLen, alignMatrix), vector<int>({ 2, 3, 5, 8, 24 })); */
}

TEST(GenotypingStrWithTwoAlleles, TypicalReads_Computed)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    AlignMatrix alignMatrix(1);
    GraphAlignment mate = decodeGraphAlignment(0, "0[6M]", &graph);
    alignMatrix.add(decodeGraphAlignment(3, "0[3M]1[1M]1[1M]2[4M]", &graph), mate);
    alignMatrix.add(decodeGraphAlignment(3, "0[3M]1[1M]1[1M]", &graph), mate);
    alignMatrix.add(decodeGraphAlignment(0, "1[1M]1[1M]1[1M]2[4M]", &graph), mate);
    alignMatrix.add(mate, decodeGraphAlignment(0, "1[1M]1[1M]1[1M]1[1M]", &graph));

    int motifLen = 1;
    int readLen = 8;
    int fragLen = 20;
    RepeatGenotype gt = genotype(AlleleCount::kTwo, motifLen, readLen, fragLen, alignMatrix);

    RepeatGenotype expectedGt(motifLen, { 2, 12 });
    expectedGt.setShortAlleleSizeInUnitsCi(2, 17);
    expectedGt.setLongAlleleSizeInUnitsCi(2, 73);
    EXPECT_EQ(gt, expectedGt);
}