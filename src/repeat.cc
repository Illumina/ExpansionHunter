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

#include "include/repeat.h"

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
#include "include/repeat_length.h"
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

void Repeat::AsPtree(ptree &allele_node) const {
  allele_node.put<string>("Size", std::to_string(size));

  if (supported_by == SupportType::kInrepeat ||
      supported_by == SupportType::kFlanking) {
    const string ci_encoding =
        std::to_string(size_ci_lower) + "," + std::to_string(size_ci_upper);
    allele_node.put<string>("CI", ci_encoding);
  }

  allele_node.put<string>("Source", readtypeToStr.at(supported_by));
  allele_node.put<size_t>("NumSupportingReads", num_supporting_reads);
}

static bool
AddConfusionCountsNode(const std::string &labelStr,
                       boost::property_tree::ptree &hunterEleNode,
                       const vector<Region> &confusionRegionTable,
                       const vector<size_t> &confusionRegionInRepeatCountVec) {
  boost::property_tree::ptree confusionCountsNode;
  size_t countInd(0);
  const bool countVecEmpty(confusionRegionInRepeatCountVec.empty());

  assert(countVecEmpty || (confusionRegionInRepeatCountVec.size() ==
                           confusionRegionTable.size()));

  for (const Region &confusionRegion : confusionRegionTable) {
    const size_t count(
        countVecEmpty ? 0 : confusionRegionInRepeatCountVec[countInd++]);
    confusionCountsNode.put<size_t>(confusionRegion.AsString(), count);
  }

  hunterEleNode.put_child(labelStr, confusionCountsNode);

  return true;
}

static bool CompareBySize(const Repeat &a1, const Repeat &a2) {
  return a1.size < a2.size;
}

void AsPtree(const Parameters &parameters, ptree &region_node,
             vector<Repeat> repeats, const RepeatSpec &region_info,
             const size_t num_irrs, const size_t num_unaligned_irrs,
             const size_t num_anchored_irrs,
             const vector<size_t> &off_target_irr_counts, vector<int> &genotype,
             const vector<array<int, 3>> &genotype_support) {
  region_node.put<string>("RepeatId", region_info.repeat_id);
  const string unit_encoding = boost::algorithm::join(region_info.units, "/");
  region_node.put<string>("RepeatUnit", unit_encoding);
  region_node.put<string>("TargetRegion", region_info.target_region.AsString());
  vector<string> genotype_encoding_vec;
  vector<string> genotype_ci_encoding_vec;

  vector<Repeat> genotype_repeats;

  for (int size : genotype) {
    genotype_encoding_vec.push_back(std::to_string(size));
    bool repeat_found = false;
    for (const Repeat &repeat : repeats) {
      if (repeat.size == size) {
        genotype_repeats.push_back(repeat);
        string ci = ".";
        if (repeat.supported_by == Repeat::SupportType::kFlanking ||
            repeat.supported_by == Repeat::SupportType::kInrepeat) {
          ci = std::to_string(repeat.size_ci_lower) + "-" +
               std::to_string(repeat.size_ci_upper);
        }
        genotype_ci_encoding_vec.push_back(ci);
        repeat_found = true;
        break;
      }
    }
    if (!repeat_found) {
      throw std::runtime_error("ERROR: Could not find " + std::to_string(size) +
                               " among repeats of " + region_info.repeat_id);
    }
  }

  if (genotype_repeats.size() == 2 &&
      genotype_repeats[0].supported_by == genotype_repeats[1].supported_by) {
    const Repeat &repeat = genotype_repeats[0];

    if (repeat.supported_by == Repeat::SupportType::kInrepeat) {

      const size_t unit_len = region_info.units[0].length();
      const double haplotype_depth = parameters.depth() / 2;

      // Calculate CI for the short allele.
      size_t short_allele_size, short_allele_size_ci_lower,
          short_allele_size_ci_upper;

      EstimateRepeatLen(num_irrs / 2, parameters.read_len(), haplotype_depth,
                        short_allele_size, short_allele_size_ci_lower,
                        short_allele_size_ci_upper);

      short_allele_size /= unit_len;
      short_allele_size_ci_lower /= unit_len;
      short_allele_size_ci_upper /= unit_len;

      // Calculate CI for the long allele.
      size_t long_allele_size, long_allele_size_ci_lower,
          long_allele_size_ci_upper;

      EstimateRepeatLen(num_irrs, parameters.read_len(), haplotype_depth,
                        long_allele_size, long_allele_size_ci_lower,
                        long_allele_size_ci_upper);

      long_allele_size /= unit_len;
      long_allele_size_ci_lower /= unit_len;
      long_allele_size_ci_upper /= unit_len;

      genotype_encoding_vec = {std::to_string(short_allele_size),
                               std::to_string(long_allele_size)};

      const string short_allele_size_ci =
          std::to_string(parameters.read_len()) + "-" +
          std::to_string(short_allele_size_ci_upper);
      const string long_allele_size_ci =
          std::to_string(short_allele_size_ci_lower) + "-" +
          std::to_string(long_allele_size_ci_upper);
      genotype_ci_encoding_vec = {short_allele_size_ci, long_allele_size_ci};

    } else if (repeat.supported_by == Repeat::SupportType::kFlanking) {
      // If both alleles are flanking, they are assumed to have the same
      // lengths and CIs.
    }
  }

  const string genotype_encoding =
      boost::algorithm::join(genotype_encoding_vec, ",");
  region_node.put<string>("Genotype", genotype_encoding);
  const string genotype_ci_encoding =
      boost::algorithm::join(genotype_ci_encoding_vec, ",");
  region_node.put<string>("GenotypeCi", genotype_ci_encoding);

  cerr << genotype_encoding << "\t" << genotype_ci_encoding << endl;

  vector<string> haplotype_support_encodings;
  for (const auto &haplotype_support : genotype_support) {
    std::stringstream stream;
    std::copy(haplotype_support.begin(), haplotype_support.end(),
              std::ostream_iterator<int>(stream, "-"));
    string haplotype_support_encoding = stream.str();
    haplotype_support_encoding = haplotype_support_encoding.substr(
        0, haplotype_support_encoding.length() - 1);
    haplotype_support_encodings.push_back(haplotype_support_encoding);
  }
  const string genotype_support_encoding =
      boost::algorithm::join(haplotype_support_encodings, ",");
  region_node.put<string>("GenotypeSupport", genotype_support_encoding);
  region_node.put<size_t>("AnchoredIrrCount", num_anchored_irrs);

  AddConfusionCountsNode("OffTargetRegionIrrCounts", region_node,
                         region_info.offtarget_regions, off_target_irr_counts);

  region_node.put<size_t>("UnalignedIrrCount", num_unaligned_irrs);
  region_node.put<size_t>("IrrCount", num_irrs);

  size_t num_allele = 1;
  ptree repeat_sizes_node;

  std::sort(repeats.begin(), repeats.end(), CompareBySize);
  for (const Repeat &repeat : repeats) {
    const string name = "Allele" + std::to_string(num_allele);
    ptree allele_node;
    repeat.AsPtree(allele_node);
    repeat_sizes_node.put_child(name, allele_node);
    ++num_allele;
  }

  region_node.put_child("RepeatSizes", repeat_sizes_node);
}

void DumpVcf(const Parameters &options,
             const map<string, RepeatSpec> repeat_specs, const ptree &root_node,
             Outputs &outputs) {
  std::stringstream vcf_header, vcf_body;

  // clang-format off
  vcf_header
      << "##fileformat=VCFv4.1\n"
         "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
         "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n"
         "##INFO=<ID=REF,Number=1,Type=Integer,Description=\"Reference copy number\">\n"
         "##INFO=<ID=RL,Number=1,Type=Integer,Description=\"Reference length in bp\">\n"
         "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat unit in the reference orientation\">\n"
         "##INFO=<ID=REPID,Number=1,Type=String,Description=\"Repeat identifier from the input specification file\">\n"
         "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
         "##FORMAT=<ID=SO,Number=1,Type=String,Description=\"Type of reads that support the allele; can be SPANNING, FLANKING, or INREPEAT meaning that the reads span, flank, or are fully contained in the repeat\">\n"
         "##FORMAT=<ID=CN,Number=1,Type=String,Description=\"Allele copy number\">\n"
         "##FORMAT=<ID=CI,Number=1,Type=String,Description=\"Confidence interval for CN\">\n"
         "##FORMAT=<ID=AD_FL,Number=1,Type=String,Description=\"Number of flanking reads consistent with the allele\">\n"
         "##FORMAT=<ID=AD_SP,Number=1,Type=String,Description=\"Number of spanning reads consistent with the allele\">\n"
         "##FORMAT=<ID=AD_IR,Number=1,Type=String,Description=\"Number of in-repeat reads consistent with the allele\">\n";
  // clang-format on

  std::set<size_t> alt_sizes;
  for (const ptree::value_type &name_region : root_node) {
    const string region_id = name_region.first;
    const ptree &region_node = name_region.second;

    if (region_id == "BamStats") {
      continue;
    }

    const string region_encoding = region_node.get<string>("TargetRegion");
    const Region region(region_encoding);
    const RepeatSpec &region_info = repeat_specs.at(region_encoding);
    const string ref(1, region_info.LeftFlankBase());
    const size_t unit_len = region_info.units[0].length();
    const size_t reference_size = region_info.ref_seq.length() / unit_len;

    const ptree &alleles_node = region_node.get_child("RepeatSizes");
    const string motif = boost::algorithm::join(region_info.units, "/");

    string alt;
    size_t genotype_num = 0;

    const string genotype_encoding = region_node.get<string>("Genotype");
    vector<string> genotype;
    boost::split(genotype, genotype_encoding, boost::is_any_of(","));

    const string genotype_ci_encoding = region_node.get<string>("GenotypeCi");
    vector<string> genotype_ci;
    boost::split(genotype_ci, genotype_ci_encoding, boost::is_any_of(","));

    const string genotype_support_encoding =
        region_node.get<string>("GenotypeSupport");
    vector<string> genotype_support;
    boost::split(genotype_support, genotype_support_encoding,
                 boost::is_any_of(","));

    if (genotype.size() != genotype_ci.size() ||
        genotype_ci.size() != genotype_support.size()) {
      throw std::runtime_error("ERROR: Inconsistent number of elements in "
                               "Genotype, GenotypeCi, and GenotypeSupport");
    }

    string format_gt, format_so, format_cn, format_ci, format_ad_sp,
        format_ad_fl, format_ad_ir;

    for (int i = 0; i != genotype.size(); ++i) {
      const int allele_size = std::stoi(genotype[i]);
      const string allele_ci = genotype_ci[i];
      const string allele_support_enc = genotype_support[i];

      vector<string> allele_support;
      boost::split(allele_support, allele_support_enc, boost::is_any_of("-"));
      assert(allele_support.size() == 3);

      ptree const *repeat_node = nullptr;

      // Only homozygous in-repeat and flanking alleles would two "-" in CI
      // encoding.
      long num_dashes = std::count(genotype_ci_encoding.begin(),
                                   genotype_ci_encoding.end(), '-');
      if (num_dashes == 2) { // homozygous in-repeat or flanking repeat.
        for (const ptree::value_type &name_repeat : alleles_node) {
          const string source = name_repeat.second.get<string>("Source");
          if (source == "INREPEAT" || source == "FLANKING") {
            repeat_node = &name_repeat.second;
            break;
          }
        }
      } else {
        for (const ptree::value_type &name_repeat : alleles_node) {
          const size_t repeat_size = name_repeat.second.get<size_t>("Size");
          if (allele_size == repeat_size) {
            repeat_node = &name_repeat.second;
            break;
          }
        }
      }
      if (!repeat_node) {
        throw std::runtime_error("ERROR: Can't find repeat of size " +
                                 std::to_string(allele_size));
      }

      const size_t allele_len = allele_size * unit_len;

      string source = repeat_node->get<string>("Source");

      if (allele_size != reference_size) {
        alt_sizes.insert(allele_size);
        if (!alt.empty()) {
          alt += ",";
        }

        if (!format_gt.empty()) {
          format_gt += "/";
          format_so += "/";
          format_cn += "/";
          format_ci += "/";
          format_ad_sp += "/";
          format_ad_fl += "/";
          format_ad_ir += "/";
        }
        alt += "<STR" + std::to_string(allele_size) + ">";
        ++genotype_num;
        format_gt += std::to_string(genotype_num);
        format_so += source;
        format_cn += std::to_string(allele_size);
        format_ci += allele_ci;
        format_ad_sp += allele_support[0];
        format_ad_fl += allele_support[1];
        format_ad_ir += allele_support[2];
      } else {
        if (!format_gt.empty()) {
          format_gt = "/" + format_gt;
          format_so = "/" + format_so;
          format_cn = "/" + format_cn;
          format_ci = "/" + format_ci;
          format_ad_sp = "/" + format_ad_sp;
          format_ad_fl = "/" + format_ad_fl;
          format_ad_ir = "/" + format_ad_ir;
        }
        format_gt = "0" + format_gt;
        format_so = source + format_so;
        format_cn = std::to_string(allele_size) + format_cn;
        format_ci = allele_ci + format_ci;
      }
    }

    const string info = "SVTYPE=STR;END=" + std::to_string(region.end()) +
                        ";REF=" + std::to_string(reference_size) + ";RL=" +
                        std::to_string(reference_size * unit_len) + ";RU=" +
                        motif + ";REPID=" + region_id;
    if (alt.empty()) {
      alt = ".";
    }

    vcf_body << region.chrom() << "\t" << region.start() - 1 << "\t.\t" << ref
             << "\t" << alt << "\t.\tPASS\t" << info
             << "\tGT:SO:CN:CI:AD_FL:AD_SP:AD_FL:AD_IR\t" << format_gt << ":"
             << format_so << ":" << format_cn << ":" << format_ci << ":"
             << format_ad_sp << ":" << format_ad_fl << ":" << format_ad_ir
             << endl;
  }

  for (size_t size : alt_sizes) {
    vcf_header << "##ALT=<ID=STR" << size
               << ",Description=\"Allele comprised of " << size
               << " repeat units\">\n";
  }
  vcf_header << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
             << options.sample_name() << endl;
  outputs.vcf() << vcf_header.str() << vcf_body.str();
}

void CoalesceFlankingReads(const RepeatSpec &repeat_spec,
                           vector<Repeat> &repeats,
                           vector<RepeatAlign> *flanking_repaligns,
                           const size_t read_len, const double hap_depth,
                           size_t motif_len,
                           const vector<vector<string>> &units_shifts,
                           size_t min_baseq, double min_wp_score) {
  const string &left_flank = repeat_spec.left_flank;
  const string &right_flank = repeat_spec.right_flank;

  size_t longest_spanning = 0;
  for (const Repeat &repeat : repeats) {
    if (repeat.supported_by == Repeat::SupportType::kSpanning) {
      if (repeat.size > longest_spanning) {
        longest_spanning = repeat.size;
      }
    }
  }

  cerr << "\t[Longest spanning allele has size " << longest_spanning << "]"
       << endl;

  bool good_repeat_exists = false;
  size_t num_reads_from_unseen_allele = 0;
  size_t longest_flanking = 0;

  cerr << "\t[There are " << flanking_repaligns->size() << " flanking reads]"
       << endl;

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

        const size_t piece_start =
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

        const size_t piece_end =
            rep_align.right_flank_len + longest_spanning * motif_len;
        piece_bases = rep_align.read.bases.substr(
            0, rep_align.read.bases.length() - piece_end);
        piece_quals = rep_align.read.quals.substr(
            0, rep_align.read.bases.length() - piece_end);
        const size_t unit_length = units_shifts[0][0].length();
        const size_t offset =
            (unit_length - piece_bases.length() % unit_length) % unit_length;
        const vector<string> &units = units_shifts[offset];
        piece_wp_score =
            MatchRepeat(units, piece_bases, piece_quals, min_baseq);
      }

      if (0.7 > flank_wp || flank_wp > 1.0) {
        cerr << "[WARNING: flank_wp = " << flank_wp << "]" << endl;
      }

      piece_wp_score /= piece_bases.length();

      if (piece_wp_score >= min_wp_score && flank_wp >= min_wp_score) {
        good_flanking_repaligns.push_back(rep_align);
        good_repeat_exists = true;
        if (rep_align.size > longest_flanking) {
          longest_flanking = rep_align.size;
        }
      } else {
        cerr << "\t[Discarding flanking read " << rep_align.read.name << " "
             << rep_align.read.bases << "]" << endl;
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

    cerr << "\t[Found " << num_reads_from_unseen_allele
         << " flanking reads longer with long repeat]" << endl;
    cerr << "\t[longest_flanking = " << longest_flanking << "]" << endl;

    size_t len_estimate = 0;
    size_t lower_bound = 0;
    size_t upper_bound = 0;

    // Haplotype depth should be twice as high because flanking reads
    // are coming from both flanks.
    EstimateRepeatLen(num_reads_from_unseen_allele, read_len, 2 * hap_depth,
                      len_estimate, lower_bound, upper_bound);

    // estimateRepeatLen adds read_len to size estimates so
    // we need to subtract it.
    len_estimate -= read_len;
    lower_bound -= read_len;
    upper_bound -= read_len;

    len_estimate = len_estimate / motif_len + longest_spanning + 1;
    lower_bound = lower_bound / motif_len + longest_spanning + 1;
    upper_bound = upper_bound / motif_len + longest_spanning + 1;

    // Repeat must be at least at long as the longest flanking read.
    lower_bound = std::max(lower_bound, longest_flanking);
    len_estimate = std::max(len_estimate, longest_flanking);
    upper_bound = std::max(upper_bound, longest_flanking);

    // Repeat estimated from flanking reads cannot be longer than the read
    // length.
    const size_t num_rep_in_read = read_len / motif_len;
    lower_bound = std::min(lower_bound, num_rep_in_read);
    len_estimate = std::min(len_estimate, num_rep_in_read);
    upper_bound = std::min(upper_bound, num_rep_in_read);

    // Make sure that size estimates have expected properties.
    if (!(lower_bound <= len_estimate && len_estimate <= upper_bound)) {
      cerr << "\t[Warning CoalesceFlankingReads: Unexpected size estimates. "
           << "Repeat size is " << len_estimate << " (LB=" << lower_bound
           << " UB=" << upper_bound << ")]" << endl;
    }

    Repeat repeat;
    repeat.supported_by = Repeat::SupportType::kFlanking;
    repeat.size = len_estimate;
    repeat.size_ci_lower = lower_bound;
    repeat.size_ci_upper = upper_bound;
    repeat.num_supporting_reads = num_reads_from_unseen_allele;
    repeat.rep_aligns = supporting_aligns;

    repeats.push_back(repeat);
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
  for (size_t i = 0; i < top.length(); ++i) {
    plot.push_back(PlotColumn(
        top[i], add_bars && std::toupper(top[i]) == bot[i] ? '|' : ' ',
        bot[i]));
  }
}

static void PlotToStream(std::ostream &ostrm, Plot &plot) {
  // Write the rows one by one.
  for (size_t i = 0; i < plot.size(); ++i) {
    ostrm << plot[i].top;
  }
  ostrm << std::endl;
  for (size_t i = 0; i < plot.size(); ++i) {
    ostrm << plot[i].mid;
  }
  ostrm << std::endl;
  for (size_t i = 0; i < plot.size(); ++i) {
    ostrm << plot[i].bot;
  }
  ostrm << std::endl;
}

static void PlotSpanningAlign(Plot &plot, const string &read_seq,
                              const string &refPrefix, const string &refSuffix,
                              const size_t prefLen, const size_t suffLen) {
  const string ref_pref =
      refPrefix.substr(refPrefix.length() - prefLen, prefLen);
  const string ref_mid = string(read_seq.length() - suffLen - prefLen, 'R');
  const string ref_suff = refSuffix.substr(0, suffLen);

  PlotGaplessAlign(plot, read_seq, ref_pref + ref_mid + ref_suff);
}

static string LowerLowqualBases(const string &bases, const string &quals,
                                size_t lowqual_cutoff) {
  assert(bases.length() == quals.length());
  string cased_bases;
  for (size_t i = 0; i != bases.length(); ++i) {
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
                        const vector<Repeat> &repeats,
                        const vector<RepeatAlign> &flanking_repaligns,
                        ostream *out) {
  const string &left_flank = repeat_spec.left_flank;
  const string &right_flank = repeat_spec.right_flank;

  *out << repeat_spec.repeat_id << ":" << endl;

  for (const Repeat &allele : repeats) {
    *out << "  " << allele.readtypeToStr.at(allele.supported_by) << "_"
         << allele.size << ":" << endl;
    for (const RepeatAlign &rep_align : allele.rep_aligns) {
      *out << "    -\n      name: \"" << rep_align.read.name << "\"" << endl;

      if (allele.supported_by == Repeat::SupportType::kSpanning ||
          allele.supported_by == Repeat::SupportType::kFlanking) {
        *out << "      align: |" << endl;
        Plot plot;
        const string cased_based = LowerLowqualBases(
            rep_align.read.bases, rep_align.read.quals, parameters.min_baseq());
        PlotGaplessAlign(plot, "        ", "        ", false);
        PlotSpanningAlign(plot, cased_based, left_flank, right_flank,
                          rep_align.left_flank_len, rep_align.right_flank_len);
        PlotToStream(*out, plot);
      } else if (allele.supported_by == Repeat::SupportType::kInrepeat) {
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
                             vector<Repeat> *repeats,
                             vector<RepeatAlign> *flanking_repaligns) {
  const size_t unit_len = repeat_spec.units_shifts[0][0].length();
  std::sort(repeats->rbegin(), repeats->rend(), CompareBySize);
  const string &left_flank = repeat_spec.left_flank;
  const string &right_flank = repeat_spec.right_flank;
  const double kWpCutoff = 0.8;

  vector<RepeatAlign> filtered_flanking_repaligns;

  for (RepeatAlign &rep_align : *flanking_repaligns) {
    const string &bases = rep_align.read.bases;
    const string &quals = rep_align.read.quals;
    const size_t non_rep_len =
        rep_align.left_flank_len + rep_align.right_flank_len;
    assert(bases.length() >= non_rep_len);
    const size_t repeat_len = bases.length() - non_rep_len;

    bool found_align = false;

    for (Repeat &repeat : *repeats) {
      const size_t allele_len = repeat.size * unit_len;
      if (repeat_len > allele_len) {
        if (rep_align.left_flank_len) {
          assert(!rep_align.right_flank_len);
          const size_t prefix_len = rep_align.left_flank_len + allele_len;
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
            cerr << "[Reasign flanking to spanning]" << endl;
            Plot plot;
            const string cased_bases =
                LowerLowqualBases(bases, quals, parameters.min_baseq());
            PlotSpanningAlign(plot, cased_bases, left_flank, right_flank,
                              rep_align.left_flank_len, bases_suffix.length());
            PlotToStream(cerr, plot);
            cerr << endl;

            found_align = true;
            rep_align.right_flank_len = bases_suffix.length();
          }
        } else if (rep_align.right_flank_len) {
          assert(!rep_align.left_flank_len);
          const size_t suffix_len = rep_align.right_flank_len + allele_len;
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
            cerr << "[Reasign flanking to spanning]" << endl;
            Plot plot;
            const string cased_bases =
                LowerLowqualBases(bases, quals, parameters.min_baseq());
            PlotSpanningAlign(plot, cased_bases, left_flank, right_flank,
                              bases_prefix.length(), rep_align.right_flank_len);
            PlotToStream(cerr, plot);
            cerr << endl;

            found_align = true;
            rep_align.left_flank_len = bases_prefix.length();
          }
        }
        if (found_align) {
          rep_align.type = RepeatAlign::Type::kSpanning;
          rep_align.size = repeat.size;
          repeat.rep_aligns.push_back(rep_align);
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
