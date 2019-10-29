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

#include "filtering/OrientationPredictor.hh"

#include <string>

#include "gtest/gtest.h"

#include "graphalign/GappedAligner.hh"
#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphcore/Path.hh"
#include "graphutils/SequenceOperations.hh"

#include "input/RegionGraph.hh"
#include "locus_spec/LocusSpec.hh"

using graphtools::reverseComplement;
using std::string;

using namespace ehunter;

TEST(PredictingQueryOrientation, TypicalQueries_Classified)
{
    graphtools::Graph graph = graphtools::makeStrGraph("TAAT", "CCG", "CCTTA");

    OrientationPredictor orientationPredictor(&graph);

    const string read = "ATCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCTTA";

    EXPECT_EQ(OrientationPrediction::kAlignsInOriginalOrientation, orientationPredictor.predict(read));

    EXPECT_EQ(
        OrientationPrediction::kAlignsInOppositeOrientation, orientationPredictor.predict(reverseComplement(read)));

    const string homopolymer(150, 'A');
    EXPECT_EQ(OrientationPrediction::kDoesNotAlign, orientationPredictor.predict(homopolymer));
}
