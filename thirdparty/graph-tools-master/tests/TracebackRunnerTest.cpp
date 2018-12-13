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

#include "graphalign/TracebackRunner.hh"

#include "gtest/gtest.h"

using std::string;

using namespace graphtools;

TEST(PerformingTraceback, NeedlemanWunschMatrixFromCoreBaseAlignment_Traced)
{
    TracebackMatrix matrix("S/0   D/-2 D/-4 D/-6\n"
                           "I/-2  M/1  D/-1 D/-3\n"
                           "I/-4  M/-1 X/0  M/-2\n"
                           "I/-6  M/-3 I/-2 M/-1\n"
                           "I/-8  I/-5 M/-4 M/-1");

    const string query = "AAAC";
    const string reference = "AGC";

    TracebackRunner traceback_runner(matrix);
    Alignment alignment = traceback_runner.runTraceback(4, 3);

    Alignment expected_alignment(0, "1M1X1I1M");
    EXPECT_EQ(expected_alignment, alignment);
}

TEST(PerformingTraceback, LocalAlignmentOfCoreBases_Traced)
{
    // GGAT-CGAA
    //   || |
    //  CATAC
    TracebackMatrix matrix("S/0 S/0 S/0 S/0  S/0 S/0\n"
                           "S/0 S/0 S/0 S/0  S/0 S/0\n"
                           "S/0 S/0 S/0 S/0  S/0 S/0\n"
                           "S/0 S/0 M/5 D/1  M/5 D/1\n"
                           "S/0 S/0 I/1 M/10 D/6 D/2\n"
                           "S/0 M/5 D/1 I/6  M/7 M/11\n"
                           "S/0 I/1 M/2 I/2  M/3 I/7\n"
                           "S/0 S/0 M/5 D/1  M/7 I/3\n"
                           "S/0 S/0 M/5 M/2  M/2 M/4");

    const string query = "GGATCGAA";
    const string reference = "CATAC";

    TracebackRunner traceback_runner(matrix);
    Alignment alignment = traceback_runner.runTraceback(5, 5);

    Alignment expected_alignment(1, "2S2M1D1M3S");
    EXPECT_EQ(expected_alignment, alignment);
}
