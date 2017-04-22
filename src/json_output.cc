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

#include "include/json_output.h"

#include <iostream>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "third_party/json/json.hpp"
#include <boost/algorithm/string/join.hpp>

#include "common/parameters.h"
#include "common/repeat_spec.h"

using std::vector;
using std::ostream;
using std::string;
using std::map;
using std::cerr;
using std::endl;

using json = nlohmann::json;

void WriteJson(const Parameters &parameters,
               const map<string, RepeatSpec> &repeat_specs,
               const vector<RegionFindings> &sample_findings, ostream &out) {
  json results_json = json({});

  for (const auto &region_findings : sample_findings) {
    const RepeatSpec &repeat_spec = repeat_specs.at(region_findings.region_id);
    const string unit_encoding = boost::algorithm::join(repeat_spec.units, "/");

    // Encode genotype.
    vector<string> genotype_encoding_vec;
    for (int size : region_findings.genotype) {
      genotype_encoding_vec.push_back(std::to_string(size));
    }

    // Encode genotype CI.
    const string genotype_ci_encoding =
        boost::algorithm::join(region_findings.genotype_ci, "/");
    const string genotype_encoding =
        boost::algorithm::join(genotype_encoding_vec, "/");

    // Encode support vector.
    vector<string> genotype_support_encoding_vec;
    for (const auto &haplotype_support : region_findings.genotype_support) {
      genotype_support_encoding_vec.push_back(haplotype_support.ToString());
    }
    const string genotype_support_encoding =
        boost::algorithm::join(genotype_support_encoding_vec, "/");

    results_json[region_findings.region_id] = {
        {"RepeatId", repeat_spec.repeat_id},
        {"RepeatUnit", unit_encoding},
        {"TargetRegion", repeat_spec.target_region.ToString()},
        {"Genotype", genotype_encoding},
        {"GenotypeCi", genotype_ci_encoding},
        {"GenotypeSupport", genotype_support_encoding},
        {"AnchoredIrrCount", region_findings.num_anchored_irrs},
        {"UnalignedIrrCount", region_findings.num_unaligned_irrs},
        {"IrrCount", region_findings.num_irrs}};

    // Add offtarget read counts if they exist.
    if (!repeat_spec.offtarget_regions.empty()) {
      assert(repeat_spec.offtarget_regions.size() ==
             region_findings.offtarget_irr_counts.size());
      auto &offtarget_section =
          results_json[region_findings.region_id]["OffTargetRegionIrrCounts"];
      for (int i = 0; i != repeat_spec.offtarget_regions.size(); ++i) {
        offtarget_section[repeat_spec.offtarget_regions[i].ToString()] =
            region_findings.offtarget_irr_counts[i];
      }
    }

    // Add detected repeats.
    if (!region_findings.repeats.empty()) {
      auto &repeat_section =
          results_json[region_findings.region_id]["RepeatSizes"];
      int num_repeat = 1;
      vector<Repeat> repeats = region_findings.repeats;
      std::sort(repeats.begin(), repeats.end(), CompareRepeatBySize);
      for (const Repeat &repeat : repeats) {
        const string repeat_id = "Repeat" + std::to_string(num_repeat);
        repeat_section[repeat_id] = {
            {"Size", repeat.size},
            {"Source", repeat.readtypeToStr.at(repeat.supported_by)},
            {"NumSupportingReads", repeat.num_supporting_reads}};

        ++num_repeat;
      }
    }
  }

  results_json["BamStats"]["ReadLength"] = parameters.read_len();
  results_json["BamStats"]["MedianDepth"] = parameters.depth();

  out << results_json.dump(4);
}

/*

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

      const int unit_len = region_info.units[0].length();
      const double haplotype_depth = parameters.depth() / 2;

      // Calculate CI for the short allele.
      int short_allele_size, short_allele_size_ci_lower,
          short_allele_size_ci_upper;

      EstimateRepeatLen(num_irrs / 2, parameters.read_len(), haplotype_depth,
                        short_allele_size, short_allele_size_ci_lower,
                        short_allele_size_ci_upper);

      short_allele_size /= unit_len;
      short_allele_size_ci_lower /= unit_len;
      short_allele_size_ci_upper /= unit_len;

      // Calculate CI for the long allele.
      int long_allele_size, long_allele_size_ci_lower,
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

 */
