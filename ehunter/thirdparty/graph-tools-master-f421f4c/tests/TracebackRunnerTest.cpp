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
