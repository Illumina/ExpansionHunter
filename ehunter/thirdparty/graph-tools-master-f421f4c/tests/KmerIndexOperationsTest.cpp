//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
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

#include "graphalign/KmerIndexOperations.hh"

#include "gtest/gtest.h"

#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphutils/SequenceOperations.hh"

using namespace graphtools;

TEST(StrandClassification, TypicalSequence_StrandDetermined)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    KmerIndex kmer_index(graph, 3);
    EXPECT_TRUE(checkIfForwardOriented(kmer_index, "CCGCCGCCGCCG"));
    EXPECT_FALSE(checkIfForwardOriented(kmer_index, reverseComplement("CCGCCGCCGCCG")));
    EXPECT_TRUE(checkIfForwardOriented(kmer_index, "CCGACGCCTCCG"));
    EXPECT_FALSE(checkIfForwardOriented(kmer_index, reverseComplement("CCGACGCCTCCG")));
}
