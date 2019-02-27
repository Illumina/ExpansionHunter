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

#include "graphalign/LinearAlignment.hh"

#include "gtest/gtest.h"

using std::list;
using std::string;

using namespace graphtools;

TEST(AlignmentInitialization, TypicalCigarString_AlignmentCreated)
{
    Alignment alignment(3u, "3M1X2N2D2M3I1M10S");
    const list<Operation> operations = { Operation("3M"), Operation("1X"), Operation("2N"), Operation("2D"),
                                         Operation("2M"), Operation("3I"), Operation("1M"), Operation("10S") };

    Alignment expected_alignment(3u, operations);
    ASSERT_EQ(expected_alignment, alignment);
}

TEST(GettingAlignmentSpans, TypicalAlignment_QueryAndreferenceSpansObtained)
{
    Alignment alignment(3u, "3M1X2M2D2M3I1M10S");
    EXPECT_EQ(22u, alignment.queryLength());
    EXPECT_EQ(11u, alignment.referenceLength());
}

TEST(EncodingAlignment, TypicalAlignment_CigarStringObtained)
{
    const string cigar_string = "3M1X2N2D2M3I1M10S";
    Alignment alignment(3u, cigar_string);

    EXPECT_EQ(cigar_string, alignment.generateCigar());
}

TEST(SplittingAlignment, SplitPositionBetweenOperations_PefixAndSuffixAlignments)
{
    // query: -AATTCGTT--TTGGGTCCCCCCCCCC
    //           ||| ||  ||   |
    //   ref: CCCTTCCNNAATT---T----------

    const string cigar_string = "2S3M1X2N2D2M3I1M10S";
    Alignment alignment(3, cigar_string);

    Alignment suffix = alignment.splitAtReferencePosition(13);

    Alignment expected_prefix(3, "2S3M1X2N2D2M3I");
    Alignment expected_suffix(13, "1M10S");
    EXPECT_EQ(expected_prefix, alignment);
    EXPECT_EQ(expected_suffix, suffix);
}

TEST(SplittingAlignment, OperationOverlapsSplitPosition_PrefixAndSuffixAlignments)
{
    Alignment alignment(0, "4M1I4M");

    Alignment suffix = alignment.splitAtReferencePosition(5);

    Alignment expected_prefix(0, "4M1I1M");
    Alignment expected_suffix(5, "3M");
    EXPECT_EQ(alignment, expected_prefix);
    EXPECT_EQ(suffix, expected_suffix);
}

TEST(SplittingAlignment, TypicalAlignments_AlignmentStatsUpdated)
{
    // query: -AATTCGTT--T TGGGTCCCCCCCCCC
    //           ||| ||  | |   |
    //   ref: CCCTTCCNNAAT T---T----------

    Alignment alignment(3, "2S3M1X2M2D2M3I1M10S");
    alignment.splitAtReferencePosition(12);

    EXPECT_EQ(6u, alignment.numMatched());
    EXPECT_EQ(1u, alignment.numMismatched());
    EXPECT_EQ(2u, alignment.numClipped());
    EXPECT_EQ(0u, alignment.numInserted());
    EXPECT_EQ(2u, alignment.numDeleted());
}

TEST(SplittingAlignment, InvalidSplitPosition_ExceptionThrown)
{
    Alignment alignment(0, "3M");
    EXPECT_ANY_THROW(alignment.splitAtReferencePosition(0));
    EXPECT_ANY_THROW(alignment.splitAtReferencePosition(3));
    EXPECT_ANY_THROW(alignment.splitAtReferencePosition(4));
}

TEST(ReversingAlignment, TypicalAlignment_ReversedAlignment)
{
    //   AAC-TCGA
    //     |  ||
    // TTTTCG-CGCC
    Alignment alignment(4, "2S1M1D1I2M1S");

    const size_t reference_length = 10;
    alignment.reverse(reference_length);

    //  AGCT-CAA
    //   ||  |
    // CCGC-GCTTTT
    Alignment expected_alignment(2, "1S2M1I1D1M2S");
    EXPECT_EQ(expected_alignment, alignment);
}
