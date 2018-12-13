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

#include "graphalign/PinnedAligner.hh"

#include "gtest/gtest.h"

using std::string;

using namespace graphtools;

TEST(PopulatingTracebackMatrix, CoreBases_TracebackMatrixPopulated)
{
    const int32_t match_score = 1;
    const int32_t mismatch_score = -1;
    const int32_t gap_penalty = -2;
    PinnedAligner aligner(match_score, mismatch_score, gap_penalty);

    const string query = "AAAC";
    const string reference = "AGC";
    TracebackMatrix matrix = aligner.populateTracebackMatrix(reference, query);

    TracebackMatrix expected_matrix("S/0   D/-2 D/-4 D/-6\n"
                                    "I/-2  M/1  D/-1 D/-3\n"
                                    "I/-4  M/-1 X/0  X/-2\n"
                                    "I/-6  M/-3 X/-2 X/-1\n"
                                    "I/-8  I/-5 X/-4 M/-1");
    EXPECT_EQ(expected_matrix, matrix);
}

TEST(PerformingPrefixAlignment, CoreBases_Aligned)
{
    // query:     TAACTTTTGGG
    //            |  |||||
    // reference: TG-CTTTTAA

    const string query = "TAACTTTTGGG";
    const string reference = "TGCTTTTAA";

    PinnedAligner aligner(1, -1, -2);
    Alignment alignment = aligner.prefixAlign(reference, query);

    Alignment expected_alignment(0, "1M1I1X5M3S");
    EXPECT_EQ(expected_alignment, alignment);
}

TEST(PerformingPrefixAlignment, NoBasesAlign_SoftclipAlignment)
{
    const string query = "AAAAA";
    const string reference = "TGCTTTT";

    PinnedAligner aligner(1, -1, -2);
    Alignment alignment = aligner.prefixAlign(reference, query);

    Alignment expected_alignment(0, "5S");
    EXPECT_EQ(expected_alignment, alignment);
}

TEST(PerformingSuffixAlignment, CoreBases_Aligned)
{
    // TCACG-GAGA
    //   ||| |||
    //  TACGAGAG-

    const string query = "TCACGGAGA";
    const string reference = "TACGAGAG";

    PinnedAligner aligner(5, -4, -8);
    Alignment alignment = aligner.suffixAlign(reference, query);

    Alignment expected_alignment(1, "2S3M1D3M1I");
    EXPECT_EQ(expected_alignment, alignment);
}

TEST(PerformingSuffixAlignment, NoBasesAlign_SoftclipAlignment)
{
    const string query = "CGCGCG";
    const string reference = "TATATATA";

    PinnedAligner aligner(5, -4, -8);
    Alignment alignment = aligner.suffixAlign(reference, query);

    Alignment expected_alignment(8, "6S");
    EXPECT_EQ(expected_alignment, alignment);
}
