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

#include <boost/algorithm/string/join.hpp>
#include "third_party/json/json.hpp"

#include "common/parameters.h"
#include "common/repeat_spec.h"

using std::vector;
using std::ostream;
using std::string;
using std::map;
using std::cerr;
using std::endl;

using json = nlohmann::json;

namespace ehunter {
void WriteJson(const Parameters &parameters,
               const map<string, RepeatSpec> &repeat_specs,
               const vector<RegionFindings> &sample_findings, ostream &out) {
  json results_json = json({});

  for (const auto &region_findings : sample_findings) {
    const RepeatSpec &repeat_spec = repeat_specs.at(region_findings.region_id);
    const string unit_encoding = boost::algorithm::join(repeat_spec.units, "/");

    // Encode GenotypeRepeat.
    vector<string> genotype_encoding_vec, genotype_ci_encoding_vec;
    for (const RepeatAllele allele : region_findings.genotype) {
      genotype_encoding_vec.push_back(std::to_string(allele.size_));
      genotype_ci_encoding_vec.push_back(allele.ci_.ToString());
    }

    // Encode GenotypeRepeat CI.
    const string genotype_ci_encoding =
        boost::algorithm::join(genotype_ci_encoding_vec, "/");
    const string genotype_encoding =
        boost::algorithm::join(genotype_encoding_vec, "/");

    // Encode support vector.
    vector<string> genotype_support_encoding_vec;
    for (const RepeatAllele &allele : region_findings.genotype) {
      genotype_support_encoding_vec.push_back(allele.support_.ToString());
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
      for (int i = 0; i != (int)repeat_spec.offtarget_regions.size(); ++i) {
        offtarget_section[repeat_spec.offtarget_regions[i].ToString()] =
            region_findings.offtarget_irr_counts[i];
      }
    }

    // Add detected repeats.
    if (!region_findings.read_groups.empty()) {
      auto &repeat_section =
          results_json[region_findings.region_id]["RepeatSizes"];
      int num_repeat = 1;
      vector<RepeatReadGroup> read_groups = region_findings.read_groups;
      std::sort(read_groups.begin(), read_groups.end(),
                CompareReadGroupsBySize);
      for (const RepeatReadGroup &read_group : read_groups) {
        const string repeat_id = "Repeat" + std::to_string(num_repeat);
        repeat_section[repeat_id] = {
            {"Size", read_group.size},
            {"Source", kReadTypeToString.at(read_group.read_type)},
            {"NumSupportingReads", read_group.num_supporting_reads}};

        ++num_repeat;
      }
    }
  }

  results_json["BamStats"]["ReadLength"] = parameters.read_len();
  results_json["BamStats"]["MedianDepth"] = parameters.depth();

  out << results_json.dump(4);
}
}  // namespace ehunter
