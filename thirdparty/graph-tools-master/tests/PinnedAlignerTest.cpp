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
