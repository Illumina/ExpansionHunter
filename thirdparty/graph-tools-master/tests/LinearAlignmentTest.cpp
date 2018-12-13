//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
