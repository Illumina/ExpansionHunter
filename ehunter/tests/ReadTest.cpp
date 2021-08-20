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

#include "core/Read.hh"

#include "gtest/gtest.h"

using namespace ehunter;

TEST(ReadInitialization, TypicalCoreInfo_CoreInfoAddedToRead)
{
    ReadId readId("frag1", MateNumber::kSecondMate);
    Read read(readId, "ATTC", true);
    EXPECT_EQ("frag1", read.fragmentId());
    EXPECT_EQ(true, read.isReversed());
}

TEST(ReadReverseComplement, SequenceReversed_ReversedFlagReversed)
{
    ReadId readId("frag1", MateNumber::kSecondMate);
    Read read(readId, "ATTCCG", true);
    EXPECT_EQ("ATTCCG", read.sequence());
    ASSERT_TRUE(read.isReversed());
    read.reverseComplement();
    EXPECT_EQ("CGGAAT", read.sequence());
    ASSERT_FALSE(read.isReversed());
}
