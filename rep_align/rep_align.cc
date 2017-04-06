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
using std::string;
#include <vector>
using std::vector;
#include <cassert>
#include <iostream>
using std::cerr;
using std::endl;

#include "include/allele.h"
#include "purity/purity.h"
#include "rep_align/rep_align.h"

size_t CountUnitsAtOffset(const vector<string> &units, const string &bases,
                          size_t offset) {
  const size_t unit_len = units[0].length();
  size_t num_units = 0;

  while (offset + unit_len <= bases.length()) {
    for (const string &unit : units) {
      if (std::equal(bases.begin() + offset, bases.begin() + offset + unit_len,
                     unit.begin())) {
        ++num_units;
        break;
      }
    }
    offset += unit_len;
  }

  return num_units;
}

size_t GetOffsetMostUnits(const vector<string> &units, const string &bases,
                          size_t *max_unit_count) {
  const size_t unit_len = units[0].length();
  size_t offset_with_most_occurances = 0;
  *max_unit_count = 0;

  for (size_t offset = 0; offset != unit_len; ++offset) {
    size_t unit_count_at_offset = CountUnitsAtOffset(units, bases, offset);
    if (unit_count_at_offset > *max_unit_count) {
      *max_unit_count = unit_count_at_offset;
      offset_with_most_occurances = offset;
    }
  }
  return offset_with_most_occurances;
}

static const string RevComp(const string &bases) {
  string bases_rc = bases;
  string::reverse_iterator rev_iter = bases_rc.rbegin();
  char complement_base(' ');

  for (const char base : bases) {
    switch (base) {
    case 'A':
      complement_base = 'T';
      break;
    case 'C':
      complement_base = 'G';
      break;
    case 'G':
      complement_base = 'C';
      break;
    case 'T':
      complement_base = 'A';
      break;
    case '/':
      complement_base = '/';
      break;
    default:
      complement_base = 'N';
    }

    *rev_iter++ = complement_base;
  }

  return bases_rc;
}

static bool PerfectMatch(const string &target_kmer,
                         const vector<string> &units) {
  assert(target_kmer.length() == units[0].length());

  for (const string &unit : units) {
    if (unit == target_kmer) {
      return true;
    }
  }
  return false;
}

bool AlignLeftFlank(const vector<string> &units, const string &left_flank,
                    const string &bases, const string &quals,
                    size_t offset_most_units, size_t min_baseq,
                    double min_wp_score, size_t *left_flank_len,
                    double *left_flank_score) {
  const size_t unit_len = units[0].length();

  *left_flank_len = 0;
  for (size_t offset = offset_most_units; offset + unit_len < bases.length();
       offset += unit_len) {
    const string cur_kmer = bases.substr(offset, unit_len);
    if (PerfectMatch(cur_kmer, units)) {
      const string bases_pref = bases.substr(0, offset);
      const string quals_pref = quals.substr(0, offset);
      const string left_flank_pref =
          left_flank.substr(left_flank.length() - offset, offset);
      const vector<string> left_flank_pref_units = {left_flank_pref};
      *left_flank_score = MatchUnits(left_flank_pref_units, bases_pref.begin(),
                                     bases_pref.end(), quals_pref.begin(),
                                     quals_pref.end(), min_baseq);
      if (*left_flank_score / bases_pref.length() >= min_wp_score) {
        const string bases_pref_rc = RevComp(bases_pref);
        const string quals_pref_rc(quals_pref.rbegin(), quals_pref.rend());

        vector<string> units_rc;
        for (const string &unit : units) {
          units_rc.push_back(RevComp(unit));
        }
        float prefix_repeat_score =
            MatchRepeat(units_rc, bases_pref_rc, quals_pref_rc, min_baseq);

        if (*left_flank_score >= 2 + prefix_repeat_score) {
          *left_flank_len = offset;
          return true;
        }
      }
    }
  }
  *left_flank_score = 0;
  return false;
}

bool AlignRightFlank(const vector<string> &units, const string &right_flank,
                     const string &bases, const string &quals,
                     size_t offset_most_units, size_t min_baseq,
                     double min_wp_score, size_t *right_flank_len,
                     double *right_flank_score) {
  vector<string> units_rc;
  for (const string &unit : units) {
    units_rc.push_back(RevComp(unit));
  }
  const string left_flank = RevComp(right_flank);
  const string bases_rc = RevComp(bases);
  const string quals_re(quals.rbegin(), quals.rend());
  const size_t offset_most_units_re =
      (bases.length() - offset_most_units) % units[0].length();

  bool is_found = AlignLeftFlank(units_rc, left_flank, bases_rc, quals_re,
                                 offset_most_units_re, min_baseq, min_wp_score,
                                 right_flank_len, right_flank_score);
  return is_found;
}

bool IsSpanningOrFlankingRead(const Parameters &params,
                              const RepeatSpec &repeat_spec,
                              const string &bases, const string &quals,
                              RepeatAlign *rep_align) {
  const vector<string> &units = repeat_spec.units_shifts[0];
  size_t max_unit_count = 0;
  size_t offset_most_units = GetOffsetMostUnits(units, bases, &max_unit_count);

  double left_flank_score = 0;
  double right_flank_score = 0;

  const double kFlankMinWpScore = 0.7;

  const bool matches_left_flank =
      AlignLeftFlank(units, repeat_spec.left_flank, bases, quals,
                     offset_most_units, params.min_baseq(), kFlankMinWpScore,
                     &rep_align->left_flank_len, &left_flank_score);

  const bool matches_right_flank =
      AlignRightFlank(units, repeat_spec.right_flank, bases, quals,
                      offset_most_units, params.min_baseq(), kFlankMinWpScore,
                      &rep_align->right_flank_len, &right_flank_score);

  if (bases.length() < rep_align->left_flank_len + rep_align->right_flank_len) {
    return false;
  }

  // Read prefix matching left flank or empty if none detected.
  const string bases_prefix = bases.substr(0, rep_align->left_flank_len);
  const string quals_prefix = quals.substr(0, rep_align->left_flank_len);

  const size_t middle_len =
      bases.length() - rep_align->right_flank_len - rep_align->left_flank_len;
  const string bases_middle =
      bases.substr(rep_align->left_flank_len, middle_len);
  const string quals_middle =
      quals.substr(rep_align->left_flank_len, middle_len);

  const string bases_suffix = bases.substr(
      bases.length() - rep_align->right_flank_len, rep_align->right_flank_len);
  const string quals_suffix = quals.substr(
      bases.length() - rep_align->right_flank_len, rep_align->right_flank_len);

  double repeat_score =
      MatchRepeat(units, bases_middle, quals_middle, params.min_baseq());

  const size_t non_repeat_len =
      rep_align->left_flank_len + rep_align->right_flank_len;
  rep_align->size = (bases.length() - non_repeat_len) / units[0].length();
  rep_align->read.bases = bases;
  rep_align->read.quals = quals;

  const double read_wp =
      (left_flank_score + repeat_score + right_flank_score) / bases.length();

  if (read_wp >= params.min_wp()) {
    if (matches_left_flank && matches_right_flank) {
      rep_align->type = RepeatAlign::Type::kSpanning;
      return true;
    } else if (matches_left_flank || matches_right_flank) {
      rep_align->type = RepeatAlign::Type::kFlanking;
      return true;
    }
  }
  rep_align->size = 0;

  return false;
}

bool IsSpanningOrFlankingReadRc(const Parameters &params,
                                const RepeatSpec &repeat_spec,
                                const string &bases, const string &quals,
                                RepeatAlign *rep_align) {
  const bool forward_match = IsSpanningOrFlankingRead(
      params, repeat_spec, bases, quals, &(*rep_align));
  if (forward_match) {
    return true;
  }

  const string bases_rc = RevComp(bases);
  const string quals_rc(quals.rbegin(), quals.rend());

  const bool reverse_match = IsSpanningOrFlankingRead(
      params, repeat_spec, bases_rc, quals_rc, &(*rep_align));

  return reverse_match;
}

static float ScoreSpanningAlign(size_t min_baseq, double min_wp,
                                const vector<string> &units,
                                const string &left_flank,
                                const string &right_flank, const string &bases,
                                const string &quals, size_t left_flank_len,
                                size_t right_flank_len) {
  const double kFlankMinWpScore = 0.7;
  const string bases_prefix = bases.substr(0, left_flank_len);
  const string quals_prefix = quals.substr(0, right_flank_len);

  assert(bases.length() >= left_flank_len + right_flank_len);
  const size_t repeat_len = bases.length() - left_flank_len - right_flank_len;

  const string bases_repeat = bases.substr(left_flank_len, repeat_len);
  const string quals_repeat = quals.substr(left_flank_len, repeat_len);

  const string bases_suffix =
      bases.substr(bases.length() - right_flank_len, right_flank_len);
  const string quals_suffix =
      quals.substr(quals.length() - right_flank_len, right_flank_len);

  float repeat_score =
      MatchRepeat(units, bases_repeat, quals_repeat, min_baseq);

  const string left_flank_pref =
      left_flank.substr(left_flank.length() - left_flank_len, left_flank_len);
  const vector<string> left_flank_pref_units = {left_flank_pref};
  float left_flank_score = MatchUnits(
      left_flank_pref_units, bases_prefix.begin(), bases_prefix.end(),
      quals_prefix.begin(), quals_prefix.end(), min_baseq);

  const string right_flank_pref = right_flank.substr(0, right_flank_len);
  const vector<string> right_flank_pref_units = {right_flank_pref};
  float right_flank_score = MatchUnits(
      right_flank_pref_units, bases_suffix.begin(), bases_suffix.end(),
      quals_suffix.begin(), quals_suffix.end(), min_baseq);

  return (left_flank_score + repeat_score + right_flank_score) / bases.length();
}

static size_t FindTopRightFlankLen(size_t min_baseq, double min_wp,
                                   const vector<string> &units,
                                   const string &left_flank,
                                   const string &right_flank,
                                   const string &bases, const string &quals,
                                   size_t cur_size, size_t cur_left_len) {
  size_t unit_len = units[0].length();
  double top_wp = 0;
  size_t top_right_len = 0;
  size_t top_size = 0;

  for (size_t test_size = 1; test_size != cur_size + 1; ++test_size) {
    const size_t test_repeat_len = test_size * unit_len;

    assert(bases.length() >= cur_left_len + test_repeat_len);
    size_t test_right_len = bases.length() - cur_left_len - test_repeat_len;
    const float test_wp =
        ScoreSpanningAlign(min_baseq, min_wp, units, left_flank, right_flank,
                           bases, quals, cur_left_len, test_right_len);
    if (test_wp > top_wp) {
      top_size = test_size;
      top_right_len = test_right_len;
      top_wp = test_wp;
    }
  }

  return top_right_len;
}

static size_t FindTopLeftFlankLen(size_t min_baseq, double min_wp,
                                  const vector<string> &units,
                                  const string &left_flank,
                                  const string &right_flank,
                                  const string &bases, const string &quals,
                                  size_t cur_size, size_t cur_right_len) {
  size_t unit_len = units[0].length();
  double top_wp = 0;
  size_t top_left_len = 0;
  size_t top_size = 0;

  for (size_t test_size = 1; test_size != cur_size + 1; ++test_size) {
    const size_t test_repeat_len = test_size * unit_len;
    assert(bases.length() >= cur_right_len + test_repeat_len);
    size_t test_left_len = bases.length() - cur_right_len - test_repeat_len;
    const float test_wp =
        ScoreSpanningAlign(min_baseq, min_wp, units, left_flank, right_flank,
                           bases, quals, test_left_len, cur_right_len);
    if (test_wp > top_wp) {
      top_size = test_size;
      top_left_len = test_left_len;
      top_wp = test_wp;
    }
  }

  return top_left_len;
}

bool AlignRead(const Parameters &params, const RepeatSpec &repeat_spec,
               const string &bases, const string &quals,
               RepeatAlign *rep_align) {
  const bool aligns = IsSpanningOrFlankingReadRc(params, repeat_spec, bases,
                                                 quals, &(*rep_align));

  if (!aligns || rep_align->type != RepeatAlign::Type::kSpanning) {
    return aligns;
  }

  const string &oriented_bases = rep_align->read.bases;
  const string &oriented_quals = rep_align->read.quals;

  // Search for a better alignment.
  const vector<string> &units = repeat_spec.units_shifts[0];
  const size_t unit_len = units[0].length();

  const size_t top_left_len = FindTopLeftFlankLen(
      params.min_baseq(), params.min_wp(), units, repeat_spec.left_flank,
      repeat_spec.right_flank, oriented_bases, oriented_quals, rep_align->size,
      rep_align->right_flank_len);

  assert(oriented_bases.length() >= top_left_len + rep_align->right_flank_len);
  size_t cur_size =
      (oriented_bases.length() - top_left_len - rep_align->right_flank_len) /
      unit_len;

  const size_t top_right_len = FindTopRightFlankLen(
      params.min_baseq(), params.min_wp(), units, repeat_spec.left_flank,
      repeat_spec.right_flank, oriented_bases, oriented_quals, cur_size,
      top_left_len);

  cur_size =
      (oriented_bases.length() - top_left_len - top_right_len) / unit_len;

  if (top_left_len != rep_align->left_flank_len ||
      top_right_len != rep_align->right_flank_len) {
    rep_align->left_flank_len = top_left_len;
    rep_align->right_flank_len = top_right_len;
    rep_align->size = cur_size;
  }
  return true;
}
