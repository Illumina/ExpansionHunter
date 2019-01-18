//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
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

#include "alignment/SoftclippingAligner.hh"

#include <string>

#include "gtest/gtest.h"

#include "graphalign/GappedAligner.hh"
#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"

#include "input/GraphBlueprint.hh"
#include "input/RegionGraph.hh"
#include "region_spec/LocusSpecification.hh"

using graphtools::decodeGraphAlignment;
using graphtools::GappedGraphAligner;
using graphtools::Graph;
using graphtools::GraphAlignment;
using std::list;
using std::string;

using namespace ehunter;

/*
class AligningReads : public ::testing::TestWithParam<std::string>
{
};

// TODO: Throw an error if there are no valid extensions?
TEST_P(AligningReads, ReadFlankedByLowQualityBases_Aligned)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATCGATCG", "(CAG)CAACAG(CCG)", "GCTAGCTA"));
    GraphAlignmentHeuristicsParameters alignmentHeuristicsParams(14, 2, 5);
    SoftclippingAligner aligner(&graph, GetParam(), alignmentHeuristicsParams);

    const string query = "gatcgCAgCAGCAACAGCCGCCGCCGCCGgcta";
    const list<GraphAlignment> alignments = aligner.align(query);

    const list<GraphAlignment> expectedAlignments
        = { decodeGraphAlignment(6, "0[3S2M]1[3M]1[3M]2[6M]3[3M]3[3M]3[3M]3[1M6S]", &graph) };
    ASSERT_EQ(expectedAlignments, alignments);
}

TEST_P(AligningReads, RealReadEndingInManyLowQualityBases_Aligned)
{
    const string leftFlank = "AGCCCCATTCATTGCCCCGGTGCTGAGCGGCGCCGCGAGTCGGCCCGAGGCCTCCGGGGACTGCCGTGCCGGGCGGGAGACCGCCATGG"
                             "CGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTC";
    const string rightFlank = "CCTCCTCAGCTTCCTCAGCCGCCGCCGCAGGCACAGCCGCTGCTGCCTCAGCCGCAGCCGCCCCCGCCGCCGCCCCCGCCGCCACCCG"
                              "GCCCGGCTGTGGCTGAGGAGCCGCTGCACCGACCGTGAGTTTGGGCC";

    Graph graph = makeRegionGraph(decodeFeaturesFromRegex(leftFlank, "(CAG)CAACAG(CCG)", rightFlank));

    GraphAlignmentHeuristicsParameters alignmentHeuristicsParams(14, 2, 5);
    SoftclippingAligner aligner(&graph, GetParam(), alignmentHeuristicsParams);

    const string query = "CCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACaGCCGCCACCGCCGCCGCCGCCGCCGCC"
                         "GCCtCCgCAGCCtCCtCaGCCGCCGCCGCCgcCgCaGCCGCcGCcgCCgCcgcCgcc";

    const list<GraphAlignment> alignments = aligner.align(query);

    const list<GraphAlignment> expectedAlignments = { decodeGraphAlignment(
        135,
        "0[1M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]1[3M]2[6M]"
        "3[3M]3[2M1X]3[3M]3[3M]3[3M]3[3M]3[3M]3[3M]3[3M]4[5M1X4M1X17M28S]",
        &graph) };
    ASSERT_EQ(expectedAlignments, alignments);
}

INSTANTIATE_TEST_CASE_P(
    AlignerTestsInst, AligningReads, ::testing::Values(std::string("path-aligner"), std::string("dag-aligner")), );
*/