//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "gtest/gtest.h"
#include "purity/purity.h"

#include <string>
using std::string;
#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

TEST(TestUnitMatching, MatchesUnitToItself) {
  char qual_chars[] = {40, 40, 40, 40, 40, 40};
  string quals = qual_chars;
  string bases = "GGCCCC";
  vector<string> units = {"GGCCCC"};

  EXPECT_DOUBLE_EQ(
      MatchUnits(units, bases.begin(), bases.end(), quals.begin(), quals.end()),
      6.0);
}

TEST(TestUnitMatching, MatchesMultipleUnits) {
  string quals = "PPPPPP";
  string bases = "AACTCC";
  vector<string> units = {"GGCCCC", "AACTCC"};

  EXPECT_DOUBLE_EQ(
      MatchUnits(units, bases.begin(), bases.end(), quals.begin(), quals.end()),
      6.0);
}

TEST(TestUnitMatching, MatchesShortSequence) {
  string quals = "PPP";
  string bases = "AAC";
  vector<string> units = {"GGCCCC", "AACTCC"};

  EXPECT_DOUBLE_EQ(
      MatchUnits(units, bases.begin(), bases.end(), quals.begin(), quals.end()),
      3.0);
}

TEST(TestUnitMatching, MatchesLowqualBases) {
  string quals = "(PP(((";
  string bases = "AACCGG";
  vector<string> units = {"GGCCCC", "AACTCC"};

  EXPECT_DOUBLE_EQ(
      MatchUnits(units, bases.begin(), bases.end(), quals.begin(), quals.end()),
      4.5);
}

TEST(TestUnitMatching, ScoreCanBeNegative) {
  string quals = "PPPPPP";
  string bases = "AACCGG";
  vector<string> units = {"ATTTTT", "AATTTT"};

  EXPECT_DOUBLE_EQ(
      MatchUnits(units, bases.begin(), bases.end(), quals.begin(), quals.end()),
      -2.0);
}

TEST(TestRepeatMatching, RepeatMatches) {
  string quals = "PPPPPPPP";
  string bases = "ACGATGAC";
  vector<string> units = {"AAG", "ACG"};

  EXPECT_DOUBLE_EQ(MatchRepeat(units, bases, quals), 6.0);
}

TEST(TestRepeatMatching, MotifShorterByOne) {
  string quals = "PPPPPPPP";
  string bases = "ACGATGAC";
  vector<string> units = {"AAAATTT", "ACGATGA"};

  EXPECT_DOUBLE_EQ(MatchRepeat(units, bases, quals), 6.0);
}

TEST(TestRepeatMatching, EmptySequenceScoresZero) {
  string quals;
  string bases;
  vector<string> units = {"AAG", "ACG"};

  EXPECT_DOUBLE_EQ(MatchRepeat(units, bases, quals), 0);
}

TEST(TestRepeatMatching, SingletonScoresOne) {
  string quals = "B";
  string bases = "G";
  vector<string> units = {"G"};

  EXPECT_DOUBLE_EQ(MatchRepeat(units, bases, quals), 1.0);
}

TEST(TestRepeatMatching, MakeShiftedUnits) {
  vector<string> units = {"AAG", "ACG"};
  vector<vector<string>> shifted_units = {
      {"AAG", "ACG"}, {"AGA", "CGA"}, {"GAA", "GAC"}};
  ASSERT_EQ(shift_units(units), shifted_units);
}

TEST(TestRepeatMatching, RepeatMatchesWithShift) {
  vector<string> units = {"AAG", "ACG"};
  vector<vector<string>> units_shifts = shift_units(units);
  string quals = "PPPPPPPP";
  string bases = "CGACGACG";

  size_t best_offset = 0;
  ASSERT_DOUBLE_EQ(MatchRepeat(units_shifts, bases, quals, best_offset), 8.0);
}

TEST(TestRepeatMatching, CalculatesBestMatchOffset) {
  vector<string> units = {"AAG", "ACG"};
  vector<vector<string>> units_shifts = shift_units(units);
  string quals = "PPPPPPPP";
  string bases = "CGACGACG";
  size_t offset = 0;
  MatchRepeat(units_shifts, bases, quals, offset);
  ASSERT_DOUBLE_EQ(offset, 1);
}

TEST(TestRepeatMatching, RepeatMatchesReverseCompliment) {
  vector<string> units = {"AAG", "ACG"};
  vector<vector<string>> units_shifts = shift_units(units);
  string quals = "((PPPPPP";
  string bases = "AATCGTCG";
  ASSERT_DOUBLE_EQ(MatchRepeatRc(units_shifts, bases, quals), 7.0);
}
