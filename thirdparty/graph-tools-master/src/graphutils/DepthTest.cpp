//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Sai Chen <schen6@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "graphutils/DepthTest.hh"

DepthTest::DepthTest(
    int32_t expectedNumReads, double stdDeviation, double lowerSignificanceThreshold, double upperSignificanceThreshold)
    : expectedNumReads_(expectedNumReads)
    , stdDeviation_(stdDeviation)
    , lowerSignificanceThreshold_(lowerSignificanceThreshold)
    , upperSignificanceThreshold_(upperSignificanceThreshold)
    , coverageDistribution_(expectedNumReads_, stdDeviation_)
{
}

bool DepthTest::testReadCount(int32_t observedNumReads)
{
    double coverageTestPvalue = cdf(coverageDistribution_, observedNumReads);

    if (coverageTestPvalue <= 0.5)
    {
        return coverageTestPvalue < lowerSignificanceThreshold_ ? false : true;
    }
    else
    {
        coverageTestPvalue = 1 - coverageTestPvalue;
        return coverageTestPvalue < upperSignificanceThreshold_ ? false : true;
    }
}
