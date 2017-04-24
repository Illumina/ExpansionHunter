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

#include "include/vcf_output.h"

#include <iostream>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include <boost/algorithm/string/join.hpp>

using std::vector;
using std::ostream;
using std::string;
using std::map;
using std::cerr;
using std::endl;

void WriteVcf(const Parameters &parameters,
              const std::map<std::string, RepeatSpec> &repeat_specs,
              const std::vector<RegionFindings> &sample_findings,
              std::ostream &out) {

  std::stringstream vcf_header, vcf_body;
  // clang-format off
  vcf_header <<
      "##fileformat=VCFv4.1\n"
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

  std::set<int> alt_sizes;

  for (const auto &region_findings : sample_findings) {
    const string &region_id = region_findings.region_id;
    const RepeatSpec &repeat_spec = repeat_specs.at(region_id);

    const string ref_field(1, repeat_spec.LeftFlankBase());
    const int unit_len = repeat_spec.units[0].length();
    const int reference_size = repeat_spec.ref_seq.length() / unit_len;

    // const ptree &alleles_node = region_node.get_child("RepeatSizes");
    const string unit_encoding = boost::algorithm::join(repeat_spec.units, "/");

    string alt;
    int genotype_num = 0;

    if (region_findings.genotype.NumAlleles() !=
            region_findings.genotype_ci.size() ||
        region_findings.genotype_ci.size() !=
            region_findings.genotype_support.size()) {
      throw std::runtime_error("ERROR: Inconsistent number of elements in "
                               "Genotype, GenotypeCi, and GenotypeSupport");
    }

    string format_gt, format_so, format_cn, format_ci, format_ad_sp,
        format_ad_fl, format_ad_ir;

    const string genotype_ci_encoding =
        boost::algorithm::join(region_findings.genotype_ci, "/");
    const vector<int> allele_sizes =
        region_findings.genotype.ExtractAlleleSizes();
    for (int i = 0; i != allele_sizes.size(); ++i) {
      const int allele_size = allele_sizes[i];
      const string allele_ci = region_findings.genotype_ci[i];
      const AlleleSupport allele_support = region_findings.genotype_support[i];

      // Only homozygous in-repeat and flanking alleles would two "-" in CI
      // encoding.

      long num_dashes = std::count(genotype_ci_encoding.begin(),
                                   genotype_ci_encoding.end(), '-');

      RepeatReadGroup const *read_group_for_allele = nullptr;
      if (num_dashes == 2) { // homozygous in - repeat or flanking repeat.
        for (const RepeatReadGroup &repeat : region_findings.read_groups) {
          if (repeat.supported_by == RepeatReadGroup::SupportType::kInrepeat ||
              repeat.supported_by == RepeatReadGroup::SupportType::kFlanking) {
            read_group_for_allele = &repeat;
            break;
          }
        }
      } else {
        for (const RepeatReadGroup &repeat : region_findings.read_groups) {
          if (allele_size == repeat.size) {
            read_group_for_allele = &repeat;
            break;
          }
        }
      }
      if (!read_group_for_allele) {
        throw std::runtime_error("ERROR: Can't find repeat of size " +
                                 std::to_string(allele_size));
      }

      const int allele_len = allele_size * unit_len;
      const string source = read_group_for_allele->readtypeToStr.at(
          read_group_for_allele->supported_by);

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
        format_ad_sp += std::to_string(allele_support.num_spanning());
        format_ad_fl += std::to_string(allele_support.num_flanking());
        format_ad_ir += std::to_string(allele_support.num_inrepeat());
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
    const Region &region = repeat_spec.target_region;
    const string info = "SVTYPE=STR;END=" + std::to_string(region.end()) +
                        ";REF=" + std::to_string(reference_size) + ";RL=" +
                        std::to_string(reference_size * unit_len) + ";RU=" +
                        unit_encoding + ";REPID=" + region_id;
    if (alt.empty()) {
      alt = ".";
    }

    vcf_body << region.chrom() << "\t" << region.start() - 1 << "\t.\t"
             << ref_field << "\t" << alt << "\t.\tPASS\t" << info
             << "\tGT:SO:CN:CI:AD_FL:AD_SP:AD_FL:AD_IR\t" << format_gt << ":"
             << format_so << ":" << format_cn << ":" << format_ci << ":"
             << format_ad_sp << ":" << format_ad_fl << ":" << format_ad_ir
             << endl;
  }

  for (int size : alt_sizes) {
    vcf_header << "##ALT=<ID=STR" << size
               << ",Description=\"Allele comprised of " << size
               << " repeat units\">\n";
  }
  vcf_header << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
             << parameters.sample_name() << endl;
  out << vcf_header.str() << vcf_body.str();
}