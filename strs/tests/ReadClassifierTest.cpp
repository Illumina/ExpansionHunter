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

#include "strs/ReadClassifier.hh"

#include "gtest/gtest.h"

using namespace ehunter;
using std::string;

MappedRead generateRead1(int contigIndex, int64_t pos, int length)
{
    ReadId readId("frag", MateNumber::kFirstMate);
    string sequence(length, 'A');
    MappedRead read(readId, sequence, false, contigIndex, pos, 60, -1, -1, true, true, true);
    return read;
}

MappedRead generateRead2(int contigIndex, int64_t pos, int length)
{
    ReadId readId("frag", MateNumber::kSecondMate);
    string sequence(length, 'C');
    MappedRead read(readId, sequence, false, contigIndex, pos, 60, -1, -1, true, true, true);
    return read;
}

TEST(ReadClassification, TargetPair_Classified)
{
    GenomicRegion target(2, 1000, 3000);
    ReadClassifier classifier(target);

    MappedRead read = generateRead1(1, 2000, 150);
    MappedRead mate = generateRead2(1, 2300, 150);

    EXPECT_EQ(PairType::kTarget, classifier.classify(read, mate));
}