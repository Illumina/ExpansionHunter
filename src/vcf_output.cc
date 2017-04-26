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

    const string unit_encoding = boost::algorithm::join(repeat_spec.units, "/");

    string alt;
    int num_alt_allele = 0;

    string sample_gt, sample_so, sample_cn, sample_ci, sample_ad_sp,
        sample_ad_fl, sample_ad_ir;

    for (const RepeatAllele allele : region_findings.genotype) {

      const int allele_len = allele.size_ * unit_len;
      const string source = kReadTypeToString.at(allele.type_);

      if (allele.size_ != reference_size) {
        alt_sizes.insert(allele.size_);
        const bool is_hom_diploid_genotype =
            region_findings.genotype.size() == 2 &&
            region_findings.genotype[0] == region_findings.genotype[1];

        const string allele_symbol =
            "<STR" + std::to_string(allele.size_) + ">";
        if (alt.empty()) {
          alt = allele_symbol;
        } else if (!is_hom_diploid_genotype) {
          alt += "," + allele_symbol;
        }

        if (!sample_gt.empty()) {
          sample_gt += "/";
          sample_so += "/";
          sample_cn += "/";
          sample_ci += "/";
          sample_ad_sp += "/";
          sample_ad_fl += "/";
          sample_ad_ir += "/";
        }
        if (num_alt_allele != 1 || !is_hom_diploid_genotype) {
          ++num_alt_allele;
        }
        sample_gt += std::to_string(num_alt_allele);
        sample_so += source;
        sample_cn += std::to_string(allele.size_);
        sample_ci += allele.ci_.ToString();
        sample_ad_sp += std::to_string(allele.support_.num_spanning());
        sample_ad_fl += std::to_string(allele.support_.num_flanking());
        sample_ad_ir += std::to_string(allele.support_.num_inrepeat());
      } else {
        if (!sample_gt.empty()) {
          sample_gt = "/" + sample_gt;
          sample_so = "/" + sample_so;
          sample_cn = "/" + sample_cn;
          sample_ci = "/" + sample_ci;
          sample_ad_sp = "/" + sample_ad_sp;
          sample_ad_fl = "/" + sample_ad_fl;
          sample_ad_ir = "/" + sample_ad_ir;
        }
        sample_gt = "0" + sample_gt;
        sample_so = source + sample_so;
        sample_cn = std::to_string(allele.size_) + sample_cn;
        sample_ci = allele.ci_.ToString() + sample_ci;
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
             << "\tGT:SO:CN:CI:AD_SP:AD_FL:AD_IR\t" << sample_gt << ":"
             << sample_so << ":" << sample_cn << ":" << sample_ci << ":"
             << sample_ad_sp << ":" << sample_ad_fl << ":" << sample_ad_ir
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