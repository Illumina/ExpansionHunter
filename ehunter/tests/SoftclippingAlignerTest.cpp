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

#include "alignment/SoftclippingAligner.hh"

#include <string>

#include "gtest/gtest.h"

#include "graphalign/GappedAligner.hh"
#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"

#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"
#include "locus/LocusSpecification.hh"

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
