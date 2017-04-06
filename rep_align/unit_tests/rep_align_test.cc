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

#include <string>
#include <iostream>
#include <vector>

#include "common/parameters.h"
#include "common/repeat_spec.h"
#include "include/allele.h"
#include "purity/purity.h"
#include "rep_align/rep_align.h"
#include "gtest/gtest.h"

using std::cerr;
using std::endl;
using std::string;
using std::vector;

//**** CountUnitsAtOffset ****

TEST(CountUnitsAtOffset, UnitlessSeq_CountZero) {
  const vector<string> units = {"AT"};
  const string untless_seq = "CGGCGGCGGCGG";
  const size_t offset = 0;
  EXPECT_EQ(CountUnitsAtOffset(units, untless_seq, offset), 0);
}

TEST(CountUnitsAtOffset, UnitsAtOffsetZero_Counted) {
  const vector<string> units = {"AT"};
  const string seq = "ATGGATCTATA";
  const size_t offset = 0;
  EXPECT_EQ(CountUnitsAtOffset(units, seq, offset), 3);
}

TEST(CountUnitsAtOffset, MultiunitsAtOffset_Counted) {
  const vector<string> units = {"AA", "TT"};
  const string seq = "CCAACCTTCCAAC";
  const size_t offset = 2;
  EXPECT_EQ(CountUnitsAtOffset(units, seq, offset), 3);
}

TEST(CountUnitsAtOffset, UnitsAtOffset_CountedWhenOffsetSetCorrectly) {
  const vector<string> units = {"AT"};
  const string seq = "CATCCATC";
  const size_t right_offset = 1;
  const size_t wrong_offset = 0;
  EXPECT_EQ(CountUnitsAtOffset(units, seq, right_offset), 2);
  EXPECT_EQ(CountUnitsAtOffset(units, seq, wrong_offset), 0);
}

//**** GetOffsetMostUnits ****

TEST(GetOffsetMostUnits, SingleUnitString_Computed) {
  const vector<string> units = {"CGG"};
  const string seq = "AACGGAAACGGACGGAACGGAAAAA";
  size_t unit_count = 0;
  const size_t offset = GetOffsetMostUnits(units, seq, &unit_count);
  EXPECT_EQ(offset, 2);
  EXPECT_EQ(unit_count, 3);
}

TEST(GetOffsetMostUnits, MultiUnitString_Computed) {
  const vector<string> units = {"CGG", "AAA"};
  const string seq = "AACGGAAACGGACGGAACGGAAAAA";
  size_t unit_count = 0;
  const size_t offset = GetOffsetMostUnits(units, seq, &unit_count);
  EXPECT_EQ(offset, 2);
  EXPECT_EQ(unit_count, 5);
}

//**** AlignLeftFlank ****

TEST(AlignLeftFlank, PrefixHasNoUnits_Detected) {
  const string bases = "CGCGATAT";
  const string quals = "QQQQQQQQ";
  const string prefix = "CGCGCGCGCG";
  const vector<string> units = {"AT"};
  const size_t offset_most_units = 0;
  const size_t min_baseq = 20;
  const double min_wp_score = 0.9;

  double left_flank_score = 0;
  size_t left_flank_len = 0;
  ASSERT_TRUE(AlignLeftFlank(units, prefix, bases, quals, offset_most_units,
                             min_baseq, min_wp_score, &left_flank_len,
                             &left_flank_score));
  ASSERT_EQ(left_flank_len, 4);
  EXPECT_DOUBLE_EQ(left_flank_score, 4.0);
}

TEST(AlignLeftFlank, LeftFlankTooSimilarToRepeatToMatch_Rejected) {
  //                    ----RRRR
  const string bases = "ATAAATAT";
  const string quals = "QQQ(QQQQ";
  const string prefix = "AAAAATAA";
  const vector<string> units = {"AT"};
  const size_t offset_most_units = 0;
  const size_t min_baseq = 20;
  const double min_wp_score = 0.9;

  double left_flank_score = 0;
  size_t left_flank_len = 0;
  ASSERT_FALSE(AlignLeftFlank(units, prefix, bases, quals, offset_most_units,
                              min_baseq, min_wp_score, &left_flank_len,
                              &left_flank_score));
}

TEST(AlignLeftFlank, InRepeatRead_Rejected) {
  const string bases = "ATATATATATATATA";
  const string quals = "QQQQQQQQQQQQQQQ";
  const string prefix = "CGCGCGCGCGCGCGCGCGCG";
  const vector<string> units = {"AT"};
  const size_t offset_most_units = 0;
  const size_t min_baseq = 20;
  const double min_wp_score = 0.9;

  double left_flank_score = 0;
  size_t left_flank_len = 0;
  ASSERT_FALSE(AlignLeftFlank(units, prefix, bases, quals, offset_most_units,
                              min_baseq, min_wp_score, &left_flank_len,
                              &left_flank_score));
  EXPECT_EQ(left_flank_len, 0);
  EXPECT_DOUBLE_EQ(left_flank_score, 0.0);
}

//**** AlignRightFlank ****

TEST(AlignRightFlank, SuffixHasNoUnits_Detected) {
  const string bases = "ATATCGC";
  const string quals = "QQQQQQQ";
  const string prefix = "CGCGCGCGC";
  const vector<string> units = {"AT"};
  const size_t offset_most_units = 0;
  const size_t min_baseq = 20;
  const double min_wp_score = 0.9;

  double right_flank_score = 0;
  size_t repeat_offset = 0;
  ASSERT_TRUE(AlignRightFlank(units, prefix, bases, quals, offset_most_units,
                              min_baseq, min_wp_score, &repeat_offset,
                              &right_flank_score));
  ASSERT_EQ(repeat_offset, 3);
  ASSERT_DOUBLE_EQ(right_flank_score, 3.0);
}

//**** IsSpanningOrFlankingRead ****

TEST(IsSpanningOrFlankingRead, UnambigousSpanningRead_Detected) {
  //                    ------RRRRRR---
  const string bases = "CGCGCGATATATGGG";
  const string quals = "QQQQQQQQQQQQQQQ";

  RepeatSpec repeat_spec;
  repeat_spec.left_flank = "CCGCGCGCGCGCGCG";
  repeat_spec.right_flank = "GGGGGGGGGGGGGGG";
  const vector<string> units = {"AT"};
  repeat_spec.units_shifts = shift_units(units);

  Parameters params;
  params.set_min_baseq(20);
  params.set_min_wp(0.9);

  RepeatAlign ra;
  ASSERT_TRUE(IsSpanningOrFlankingRead(params, repeat_spec, bases, quals, &ra));
  ASSERT_EQ(ra.type, RepeatAlign::Type::kSpanning);
  EXPECT_EQ(ra.left_flank_len, 6);
  EXPECT_EQ(ra.right_flank_len, 3);
}

/*
TEST(IsSpanningOrFlankingRead, UnitlessRead_Rejected) {
  const string bases = "ATATATATATATATA";
  const string quals = "QQQQQQQQQQQQQQQ";
  const string left_flank = "CGCGCGCGCGGCGGCG";
  const string right_flank = "GGGGGGGGGGGGGGGGGG";
  const vector<string> units = {"CG"};
  const vector<vector<string>> units_shifts = shift_units(units);
  const size_t min_baseq = 20;
  const float min_wp_score = 0.9;

  RepeatAlign ra;
  ASSERT_FALSE(IsSpanningOrFlankingRead(units_shifts, min_baseq, min_wp_score,
                                        left_flank, right_flank, bases, quals,
                                        &ra));
  EXPECT_EQ(ra.size, 0);
}

TEST(IsSpanningOrFlankingRead, SpanningReadWrongRepeatUnit_Rejected) {
  //                    ------RRRRRR---
  const string bases = "CGCGCGATCCATGGG";
  const string quals = "QQQQQQQQQQQQQQQ";
  const string left_flank = "CCCCCCGCGCGCGCG";
  const string right_flank = "GGGGGGGGGGGGGGG";
  const vector<string> units = {"AT"};
  const vector<vector<string>> units_shifts = shift_units(units);
  const size_t min_baseq = 20;
  const float min_wp_score = 0.9;

  RepeatAlign ra;
  ASSERT_FALSE(IsSpanningOrFlankingRead(units_shifts, min_baseq, min_wp_score,
                                        left_flank, right_flank, bases, quals,
                                        &ra));
}

TEST(IsSpanningOrFlankingRead, LeftFlankingRead_Detected) {
  //                    --------RRRRRRRRRRR
  const string bases = "AAAAAAAACCGCCGCCGCC";
  const string quals = "QQQQQQQQQQQQQQQQQQQ";
  const string left_flank = "AAAAAAAAAAAAAAAAAAA";
  const string right_flank = "GGGGGGGGGGGGGGGGGGG";
  const vector<string> units = {"CCG"};
  const vector<vector<string>> units_shifts = shift_units(units);
  const size_t min_baseq = 20;
  const float min_wp_score = 0.9;

  size_t left_flank_len = 0;
  size_t right_flank_len = 0;
  RepeatAlign::Type read_type;

  RepeatAlign ra;
  ASSERT_TRUE(IsSpanningOrFlankingRead(units_shifts, min_baseq, min_wp_score,
                                       left_flank, right_flank, bases, quals,
                                       &ra));
  ASSERT_EQ(ra.type, RepeatAlign::Type::kFlanking);
  EXPECT_EQ(ra.left_flank_len, 8);
  EXPECT_EQ(ra.right_flank_len, 0);
}

TEST(IsSpanningOrFlankingRead, RightFlankingRead_Detected) {
  //                    RRRRRRRRRRRR-------
  const string bases = "GGCCCCGGCCCCGGGGGGG";
  const string quals = "QQQQQQQQQQQQQQQQQQQ";
  const string left_flank = "AAAAAAAAAAAAAAAAAAA";
  const string right_flank = "GGGGGGGGGGGGGGGGGGG";
  const vector<string> units = {"GGCCCC"};
  const vector<vector<string>> units_shifts = shift_units(units);
  const size_t min_baseq = 20;
  const float min_wp_score = 0.9;

  size_t left_flank_len = 0;
  size_t right_flank_len = 0;
  RepeatAlign::Type read_type;

  RepeatAlign ra;
  ASSERT_TRUE(IsSpanningOrFlankingRead(units_shifts, min_baseq, min_wp_score,
                                       left_flank, right_flank, bases, quals,
                                       &ra));
  ASSERT_EQ(ra.type, RepeatAlign::Type::kFlanking);
  EXPECT_EQ(ra.left_flank_len, 0);
  EXPECT_EQ(ra.right_flank_len, 7);
  EXPECT_EQ(ra.size, 2);
}

TEST(IsSpanningOrFlankingReadRc, ReverseComplimentedSpanningRead_Detected) {
  //                    ------RRRRRR---
  const string bases = "CCCATATATCGCGCG";
  const string quals = "QQQQQQQQQQQQ(((";
  const string left_flank = "CCGCGCGCGCGCGCG";
  const string right_flank = "GGGGGGGGGGGGGGG";
  const vector<string> units = {"AT"};
  const vector<vector<string>> units_shifts = shift_units(units);

  const size_t min_baseq = 20;
  const double min_wp_score = 0.9;

  RepeatAlign ra;

  ASSERT_TRUE(IsSpanningOrFlankingReadRc(units_shifts, min_baseq, min_wp_score,
                                         left_flank, right_flank, bases, quals,
                                         &ra));
  ASSERT_EQ(ra.read.bases, "CGCGCGATATATGGG");
  ASSERT_EQ(ra.read.quals, "(((QQQQQQQQQQQQ");
  ASSERT_EQ(ra.type, RepeatAlign::Type::kSpanning);
  EXPECT_EQ(ra.left_flank_len, 6);
  EXPECT_EQ(ra.right_flank_len, 3);
}

TEST(IsSpanningOrFlankingReadNew, RealCagRepeatWithRightFlankProblem_Detected) {
  //                    -----------------RRRRRRRRR-----------------
  const string bases = "AGTCCCTCAAGTCCTTCCAGCAGCAGCAACAGCCGCCGCCGCC";
  const string quals = "QQQQQQQQQ(QQQQQ(QQQQQQQQQQQQQQQQQQQQQQQQQQQ";
  const string left_flank = "CCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTC";
  const string right_flank = "CAACAGCCGCCACCGCCGCCGCCGCCGCCGCCGCCTCCTCAG";
  const vector<string> units = {"CAG"};
  const vector<vector<string>> units_shifts = shift_units(units);

  const size_t min_baseq = 20;
  const double min_wp_score = 0.9;

  RepeatAlign ra;

  ASSERT_TRUE(IsSpanningOrFlankingRead(units_shifts, min_baseq, min_wp_score,
                                       left_flank, right_flank, bases, quals,
                                       &ra));
  ASSERT_EQ(ra.type, RepeatAlign::Type::kSpanning);
  EXPECT_EQ(ra.left_flank_len, 17);
  EXPECT_EQ(ra.right_flank_len, 17);
} */