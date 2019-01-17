//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "genotyping/AllelePresenceChecker.hh"

#include <numeric>

#include "gtest/gtest.h"

#include "common/Common.hh"

using namespace ehunter;

TEST(AllelePresenceChecker, ThrowsWithIllegalParameter)
{
    ASSERT_ANY_THROW(new AllelePresenceChecker(0));
    ASSERT_ANY_THROW(new AllelePresenceChecker(15, 1));
    ASSERT_ANY_THROW(new AllelePresenceChecker(15, 0.01, -1));
    ASSERT_ANY_THROW(new AllelePresenceChecker(0));

    AllelePresenceChecker negCountGenotyper(15.0);
    ASSERT_ANY_THROW(negCountGenotyper.check(-1, 20));
}

TEST(AllelePresenceChecker, NoReads)
{
    AllelePresenceChecker checker(15.0);
    EXPECT_EQ(AllelePresenceStatus::kAbsent, checker.check(0, 0));
}

TEST(AllelePresenceChecker, AllelePresent)
{
    AllelePresenceChecker checker(15.0);
    EXPECT_EQ(AllelePresenceStatus::kPresent, checker.check(30, 30));
    EXPECT_EQ(AllelePresenceStatus::kPresent, checker.check(10, 45));
    EXPECT_EQ(AllelePresenceStatus::kPresent, checker.check(10, 0));
    EXPECT_EQ(AllelePresenceStatus::kPresent, checker.check(50, 60));
}

TEST(AllelePresenceChecker, AlleleAbsent)
{
    AllelePresenceChecker checker(15.0);
    EXPECT_EQ(AllelePresenceStatus::kAbsent, checker.check(0, 30));
    EXPECT_EQ(AllelePresenceStatus::kAbsent, checker.check(1, 60));
    EXPECT_EQ(AllelePresenceStatus::kAbsent, checker.check(1, 5));
}

TEST(AllelePresenceChecker, NoCall)
{
    AllelePresenceChecker checker(15.0);
    EXPECT_EQ(AllelePresenceStatus::kUncertain, checker.check(5, 30));
    EXPECT_EQ(AllelePresenceStatus::kUncertain, checker.check(1, 0));
}

TEST(AllelePresenceChecker, HighReads)
{
    AllelePresenceChecker checker(150.0);
    EXPECT_EQ(AllelePresenceStatus::kPresent, checker.check(100, 450));
    EXPECT_EQ(AllelePresenceStatus::kAbsent, checker.check(20, 600));
    EXPECT_EQ(AllelePresenceStatus::kUncertain, checker.check(40, 200));
}
