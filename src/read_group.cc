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

#include "include/read_group.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/property_tree/ptree.hpp>

#include "common/genomic_region.h"
#include "common/parameters.h"
#include "common/repeat_spec.h"
#include "purity/purity.h"

using std::cerr;
using std::endl;
using std::ostream;
using std::string;
using std::vector;
using std::pair;
using std::map;
using std::array;

using boost::property_tree::ptree;
// using boost::algorithm::split;
// using boost::algorithm::is_any_of;

void CoalesceFlankingReads(const RepeatSpec &repeat_spec,
                           vector<RepeatReadGroup> &repeats,
                           vector<RepeatAlign> *flanking_repaligns,
                           const int read_len, const double hap_depth,
                           int motif_len,
                           const vector<vector<string>> &units_shifts,
                           int min_baseq, double min_wp_score) {
  const string &left_flank = repeat_spec.left_flank;
  const string &right_flank = repeat_spec.right_flank;

  int longest_spanning = 0;
  for (const RepeatReadGroup &read_group : repeats) {
    if (read_group.read_type == ReadType::kSpanning) {
      if (read_group.size > longest_spanning) {
        longest_spanning = read_group.size;
      }
    }
  }

  //cerr << "\t[Longest spanning allele has size " << longest_spanning << "]"
  //     << endl;

  bool good_repeat_exists = false;
  int num_reads_from_unseen_allele = 0;
  int longest_flanking = 0;

  //cerr << "\t[There are " << flanking_repaligns->size() << " flanking reads]"
  //     << endl;

  vector<RepeatAlign> good_flanking_repaligns;

  for (const auto &rep_align : *flanking_repaligns) {
    if (rep_align.size > longest_spanning) {
      ++num_reads_from_unseen_allele;

      string piece_bases, piece_quals;

      double piece_wp_score = 0;
      double flank_wp = 0;

      if (rep_align.left_flank_len) {
        const string bases_prefix =
            rep_align.read.bases.substr(0, rep_align.left_flank_len);
        const string quals_prefix =
            rep_align.read.quals.substr(0, rep_align.left_flank_len);
        const string left_flank_pref =
            left_flank.substr(left_flank.length() - rep_align.left_flank_len,
                              rep_align.left_flank_len);
        const vector<string> left_flank_pref_units = {left_flank_pref};
        double flank_score = MatchUnits(
            left_flank_pref_units, bases_prefix.begin(), bases_prefix.end(),
            quals_prefix.begin(), quals_prefix.end(), min_baseq);
        flank_wp = flank_score / rep_align.left_flank_len;

        const int piece_start =
            rep_align.left_flank_len + longest_spanning * motif_len;
        assert(piece_start < rep_align.read.bases.length());
        piece_bases = rep_align.read.bases.substr(
            piece_start, rep_align.read.bases.length() - piece_start);
        piece_quals = rep_align.read.quals.substr(
            piece_start, rep_align.read.bases.length() - piece_start);
        const vector<string> &units = units_shifts[0];
        piece_wp_score =
            MatchRepeat(units, piece_bases, piece_quals, min_baseq);
      } else {
        assert(rep_align.right_flank_len);
        const string bases_suffix = rep_align.read.bases.substr(
            rep_align.read.bases.length() - rep_align.right_flank_len,
            rep_align.right_flank_len);
        const string quals_suffix = rep_align.read.quals.substr(
            rep_align.read.quals.length() - rep_align.right_flank_len,
            rep_align.right_flank_len);
        const string right_flank_pref =
            right_flank.substr(0, rep_align.right_flank_len);
        const vector<string> right_flank_pref_units = {right_flank_pref};
        double flank_score = MatchUnits(
            right_flank_pref_units, bases_suffix.begin(), bases_suffix.end(),
            quals_suffix.begin(), quals_suffix.end(), min_baseq);
        flank_wp = flank_score / rep_align.right_flank_len;

        const int piece_end =
            rep_align.right_flank_len + longest_spanning * motif_len;
        piece_bases = rep_align.read.bases.substr(
            0, rep_align.read.bases.length() - piece_end);
        piece_quals = rep_align.read.quals.substr(
            0, rep_align.read.bases.length() - piece_end);
        const int unit_length = units_shifts[0][0].length();
        const int offset =
            (unit_length - piece_bases.length() % unit_length) % unit_length;
        const vector<string> &units = units_shifts[offset];
        piece_wp_score =
            MatchRepeat(units, piece_bases, piece_quals, min_baseq);
      }

      //if (0.7 > flank_wp || flank_wp > 1.0) {
      //  cerr << "[WARNING: flank_wp = " << flank_wp << "]" << endl;
      //}

      piece_wp_score /= piece_bases.length();

      if (piece_wp_score >= min_wp_score && flank_wp >= min_wp_score) {
        good_flanking_repaligns.push_back(rep_align);
        good_repeat_exists = true;
        if (rep_align.size > longest_flanking) {
          longest_flanking = rep_align.size;
        }
      } else {
        //cerr << "\t[Discarding flanking read " << rep_align.read.name << " "
        //     << rep_align.read.bases << "]" << endl;
      }
    } else {
      good_flanking_repaligns.push_back(rep_align);
    }
  }

  *flanking_repaligns = good_flanking_repaligns;

  if (good_repeat_exists) {
    vector<RepeatAlign> short_aligns;
    vector<RepeatAlign> supporting_aligns;
    for (const RepeatAlign &rep_align : *flanking_repaligns) {
      if (rep_align.size > longest_spanning) {
        supporting_aligns.push_back(rep_align);
      } else {
        short_aligns.push_back(rep_align);
      }
    }
    *flanking_repaligns = short_aligns;

    //cerr << "\t[Found " << num_reads_from_unseen_allele
    //     << " flanking reads longer with long repeat]" << endl;
    //cerr << "\t[longest_flanking = " << longest_flanking << "]" << endl;

    RepeatReadGroup read_group;
    read_group.read_type = ReadType::kFlanking;
    read_group.size = longest_flanking;
    read_group.num_supporting_reads = num_reads_from_unseen_allele;
    read_group.rep_aligns = supporting_aligns;

    repeats.push_back(read_group);
  }
}

struct PlotColumn {
  PlotColumn(char t, char m, char b) {
    top = t;
    mid = m;
    bot = b;
  }
  char top;
  char mid;
  char bot;
};

typedef vector<PlotColumn> Plot;

static void PlotGaplessAlign(Plot &plot, const string &top, const string &bot,
                             const bool add_bars = true) {
  assert(top.length() == bot.length());
  for (int i = 0; i < top.length(); ++i) {
    plot.push_back(PlotColumn(
        top[i], add_bars && std::toupper(top[i]) == bot[i] ? '|' : ' ',
        bot[i]));
  }
}

static void PlotToStream(std::ostream &ostrm, Plot &plot) {
  // Write the rows one by one.
  for (int i = 0; i < plot.size(); ++i) {
    ostrm << plot[i].top;
  }
  ostrm << std::endl;
  for (int i = 0; i < plot.size(); ++i) {
    ostrm << plot[i].mid;
  }
  ostrm << std::endl;
  for (int i = 0; i < plot.size(); ++i) {
    ostrm << plot[i].bot;
  }
  ostrm << std::endl;
}

static void PlotSpanningAlign(Plot &plot, const string &read_seq,
                              const string &refPrefix, const string &refSuffix,
                              const int prefLen, const int suffLen) {
  const string ref_pref =
      refPrefix.substr(refPrefix.length() - prefLen, prefLen);
  const string ref_mid = string(read_seq.length() - suffLen - prefLen, 'R');
  const string ref_suff = refSuffix.substr(0, suffLen);

  PlotGaplessAlign(plot, read_seq, ref_pref + ref_mid + ref_suff);
}

static string LowerLowqualBases(const string &bases, const string &quals,
                                int lowqual_cutoff) {
  assert(bases.length() == quals.length());
  string cased_bases;
  for (int i = 0; i != bases.length(); ++i) {
    if (quals[i] - 33 < lowqual_cutoff) {
      cased_bases += std::tolower(bases[i]);
    } else {
      cased_bases += bases[i];
    }
  }
  return cased_bases;
}

void OutputRepeatAligns(const Parameters &parameters,
                        const RepeatSpec &repeat_spec,
                        const vector<RepeatReadGroup> &read_groups,
                        const vector<RepeatAlign> &flanking_repaligns,
                        ostream *out) {
  const string &left_flank = repeat_spec.left_flank;
  const string &right_flank = repeat_spec.right_flank;

  *out << repeat_spec.repeat_id << ":" << endl;

  for (const RepeatReadGroup &read_group : read_groups) {
    *out << "  " << kReadTypeToString.at(read_group.read_type) << "_"
         << read_group.size << ":" << endl;
    for (const RepeatAlign &rep_align : read_group.rep_aligns) {
      *out << "    -\n      name: \"" << rep_align.read.name << "\"" << endl;

      if (read_group.read_type == ReadType::kSpanning ||
          read_group.read_type == ReadType::kFlanking) {
        *out << "      align: |" << endl;
        Plot plot;
        const string cased_based = LowerLowqualBases(
            rep_align.read.bases, rep_align.read.quals, parameters.min_baseq());
        PlotGaplessAlign(plot, "        ", "        ", false);
        PlotSpanningAlign(plot, cased_based, left_flank, right_flank,
                          rep_align.left_flank_len, rep_align.right_flank_len);
        PlotToStream(*out, plot);
      } else if (read_group.read_type == ReadType::kInrepeat) {
        const string read_bases = LowerLowqualBases(
            rep_align.read.bases, rep_align.read.quals, parameters.min_baseq());
        const string mate_bases = LowerLowqualBases(
            rep_align.mate.bases, rep_align.mate.quals, parameters.min_baseq());

        if (rep_align.type == RepeatAlign::Type::kAnchored) {
          *out << "      irr: " << read_bases << endl;
          *out << "      anc: " << mate_bases << endl;
        } else if (rep_align.type == RepeatAlign::Type::kAlignedIrrPair) {
          *out << "      al_ir1: " << read_bases << endl;
          *out << "      al_ir2: " << mate_bases << endl;
        } else if (rep_align.type == RepeatAlign::Type::kUnalignedIrrPair) {
          *out << "      un_ir1: " << read_bases << endl;
          *out << "      un_ir2: " << mate_bases << endl;
        } else if (rep_align.type ==
                   RepeatAlign::Type::kUnalignedIrrSingleton) {
          *out << "      un_ir: " << read_bases << endl;
          *out << "      un_ma: " << mate_bases << endl;
        }
      } else {
        throw std::logic_error("Unknown repeat allele type");
      }
    }
  }

  if (!flanking_repaligns.empty()) {
    *out << "  FLANKING:" << endl;
    for (const RepeatAlign &rep_align : flanking_repaligns) {
      *out << "    -\n      name: \"" << rep_align.read.name << "\"" << endl;
      *out << "      align: |" << endl;
      Plot plot;
      const string cased_based = LowerLowqualBases(
          rep_align.read.bases, rep_align.read.quals, parameters.min_baseq());
      PlotGaplessAlign(plot, "        ", "        ", false);
      PlotSpanningAlign(plot, cased_based, left_flank, right_flank,
                        rep_align.left_flank_len, rep_align.right_flank_len);
      PlotToStream(*out, plot);
    }
  }
  *out << endl;
}

// Attempt to reclassify flanking reads as spanning.
void DistributeFlankingReads(const Parameters &parameters,
                             const RepeatSpec &repeat_spec,
                             vector<RepeatReadGroup> *read_groups,
                             vector<RepeatAlign> *flanking_repaligns) {
  const int unit_len = repeat_spec.units_shifts[0][0].length();
  std::sort(read_groups->rbegin(), read_groups->rend(),
            CompareReadGroupsBySize);
  const string &left_flank = repeat_spec.left_flank;
  const string &right_flank = repeat_spec.right_flank;
  const double kWpCutoff = 0.8;

  vector<RepeatAlign> filtered_flanking_repaligns;

  for (RepeatAlign &rep_align : *flanking_repaligns) {
    const string &bases = rep_align.read.bases;
    const string &quals = rep_align.read.quals;
    const int non_rep_len =
        rep_align.left_flank_len + rep_align.right_flank_len;
    assert(bases.length() >= non_rep_len);
    const int repeat_len = bases.length() - non_rep_len;

    bool found_align = false;

    for (RepeatReadGroup &read_group : *read_groups) {
      const int allele_len = read_group.size * unit_len;
      if (repeat_len > allele_len) {
        if (rep_align.left_flank_len) {
          assert(!rep_align.right_flank_len);
          const int prefix_len = rep_align.left_flank_len + allele_len;
          const string bases_suffix =
              bases.substr(prefix_len, bases.length() - prefix_len);
          const string quals_suffix =
              quals.substr(prefix_len, quals.length() - prefix_len);

          const string right_flank_ref =
              right_flank.substr(0, bases_suffix.length());
          const vector<string> right_flank_ref_units = {right_flank_ref};
          float right_flank_score = MatchUnits(
              right_flank_ref_units, bases_suffix.begin(), bases_suffix.end(),
              quals_suffix.begin(), quals_suffix.end(), parameters.min_baseq());
          if (right_flank_score / bases_suffix.length() >= kWpCutoff) {
            // cerr << "[Reasign flanking to spanning]" << endl;
            // Plot plot;
            // const string cased_bases =
            //     LowerLowqualBases(bases, quals, parameters.min_baseq());
            // PlotSpanningAlign(plot, cased_bases, left_flank, right_flank,
            //                   rep_align.left_flank_len, bases_suffix.length());
            // PlotToStream(cerr, plot);
            // cerr << endl;

            found_align = true;
            rep_align.right_flank_len = bases_suffix.length();
          }
        } else if (rep_align.right_flank_len) {
          assert(!rep_align.left_flank_len);
          const int suffix_len = rep_align.right_flank_len + allele_len;
          const string bases_prefix =
              bases.substr(0, bases.length() - suffix_len);
          const string quals_prefix =
              quals.substr(0, quals.length() - suffix_len);

          const string left_flank_ref =
              left_flank.substr(left_flank.length() - bases_prefix.length(),
                                bases_prefix.length());
          const vector<string> left_flank_ref_units = {left_flank_ref};
          double left_flank_score = MatchUnits(
              left_flank_ref_units, bases_prefix.begin(), bases_prefix.end(),
              quals_prefix.begin(), quals_prefix.end(), parameters.min_baseq());

          if (left_flank_score / bases_prefix.length() >= kWpCutoff) {
            //cerr << "[Reasign flanking to spanning]" << endl;
            //Plot plot;
            //const string cased_bases =
            //    LowerLowqualBases(bases, quals, parameters.min_baseq());
            //PlotSpanningAlign(plot, cased_bases, left_flank, right_flank,
            //                  bases_prefix.length(), rep_align.right_flank_len);
            //PlotToStream(cerr, plot);
            //cerr << endl;

            found_align = true;
            rep_align.left_flank_len = bases_prefix.length();
          }
        }
        if (found_align) {
          rep_align.type = RepeatAlign::Type::kSpanning;
          rep_align.size = read_group.size;
          read_group.rep_aligns.push_back(rep_align);
          break;
        }
      }
    }

    if (!found_align) {
      filtered_flanking_repaligns.push_back(rep_align);
    }
  }
  *flanking_repaligns = filtered_flanking_repaligns;
}
