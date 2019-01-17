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
#include "region_spec/LocusSpecification.hh"

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
        OrientationPrediction::kAlignsInReverseComplementOrientation,
        orientationPredictor.predict(reverseComplement(read)));

    const string homopolymer(150, 'A');
    EXPECT_EQ(OrientationPrediction::kDoesNotAlign, orientationPredictor.predict(homopolymer));
}
