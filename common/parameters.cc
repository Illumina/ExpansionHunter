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

#include "common/parameters.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
using boost::lexical_cast;
#include <boost/format.hpp>
using boost::format;
#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
using boost::assign::list_of;

#include <string>
using std::string;
#include <ios>
#include <iostream>
using std::cerr;
using std::endl;
#include <iomanip>
#include <stdexcept>
#include <vector>
using std::vector;

Outputs::Outputs(const string vcf_path, const string json_path,
                 const string log_path) {
  vcf_.open(vcf_path.c_str());

  if (!vcf_.is_open()) {
    throw std::runtime_error("Failed to open output file '" + vcf_path +
                             "' for writing: " + strerror(errno));
  }

  if (!json_path.empty()) {
    json_.open(json_path.c_str());

    if (!json_.is_open()) {
      throw std::runtime_error("Failed to open output JSON file '" + json_path +
                               "' for writing: " + strerror(errno));
    }
  }

  if (!log_path.empty()) {
    log_.open(log_path.c_str());

    if (!log_.is_open()) {
      throw std::runtime_error("Failed to open log file '" + log_path +
                               "' for writing: " + strerror(errno));
    }

    log_ << std::unitbuf << std::setprecision(4);
  }
}

bool CheckIndexFile(const string &bam_path) {
  vector<string> indexExtStrVec = list_of(".bai")(".csi")(".crai");

  for (const string &indexExtStr : indexExtStrVec) {
    if (boost::filesystem::exists(bam_path + indexExtStr)) {
      return true;
    }
  }

  return false;
}

void DieNotOutputFilePath(const string &output_path_str) {
  const boost::filesystem::path output_path(output_path_str);
  const boost::filesystem::path output_dir(output_path.parent_path());

  const bool is_no_dir = output_dir.empty();
  const bool is_existing_dir = boost::filesystem::is_directory(output_dir);
  const bool is_valid_fname =
      boost::filesystem::portable_posix_name(output_path.filename().string());

  if ((is_no_dir || is_existing_dir) && is_valid_fname) {
    return;
  }
  const string error_msg =
      "'" + output_path_str + "' is not a valid output path.";
  throw std::invalid_argument(error_msg);
}

bool Parameters::Load(int numArgs, char *argPtrArr[]) {
  // clang-format off
  po::options_description usage("Allowed options");
  usage.add_options()
      ("help", "Print help message")
      ("version", "Print version number")
      ("bam", po::value<string>(), "BAM file path")
      ("ref-fasta", po::value<string>(), "Reference genome file (FASTA) path")
      ("repeat-specs", po::value<string>(), "Directory containing JSON files specifying target repeat regions")
      ("vcf", po::value<string>(), "Output VCF file path")
      ("json", po::value<string>(), "Output JSON file path")
      ("log", po::value<string>(), "Output read alignment file path")
      ("region-extension-length", po::value<int>(),
          ("[Optional] How far from on/off-target regions to search for informative reads [Default " + std::to_string(region_extension_len_) + "]").c_str())
      ("min-score", po::value<float>(),
          ("[Optional] Minimum weighted matching score (0 <= x <= 1) [Default " + boost::str(format("%.3f") % min_wp_) + "]").c_str())
      ("min-baseq", po::value<int>(), ("[Optional] Minimum quality of a high-confidence base call [Default " + std::to_string(min_baseq_) + "]").c_str())
      ("min-anchor-mapq", po::value<int>(), ("[Optional] Minimum MAPQ of a read anchor [Default " + std::to_string(min_anchor_mapq_) + "]").c_str())
      ("skip-unaligned", "[Optional] Do not search for IRRs in unaligned reads")
      ("read-depth", po::value<float>(), "[Optional] Read depth")
      ("sex", po::value<string>(), "[Optional] Sex of the sample; can be either male or female [Default female]");
  // clang-format on
  po::variables_map argMap;
  po::store(po::parse_command_line(numArgs, argPtrArr, usage), argMap);
  po::notify(argMap);

  if (argMap.count("help")) {
    std::cerr << usage << std::endl;
    return false;
  }

  if (argMap.count("version")) {
    return false;
  }

  if (!argMap.count("bam")) {
    throw std::invalid_argument("--bam parameter is required");
  }

  bam_path_ = argMap["bam"].as<string>();

  if (!boost::filesystem::exists(bam_path_)) {
    throw std::invalid_argument("'" + bam_path_ + "' does not exist");
  }

  if (!CheckIndexFile(bam_path_)) {
    throw std::invalid_argument("Could not find index file for BAM '" +
                                bam_path_ + "'");
  }

  // Extract sample name.
  boost::filesystem::path boost_bam_path(bam_path_);
  sample_name_ = boost_bam_path.stem().string();

  if (!argMap.count("ref-fasta")) {
    throw std::invalid_argument("--ref-fasta parameter is required");
  }

  genome_path_ = argMap["ref-fasta"].as<string>();

  if (!boost::filesystem::exists(genome_path_)) {
    throw std::invalid_argument("'" + genome_path_ + "' does not exist");
  }

  if (!argMap.count("repeat-specs")) {
    throw std::invalid_argument("--repeat-specs parameter is required");
  }

  repeat_specs_path_ = argMap["repeat-specs"].as<string>();

  if (!boost::filesystem::exists(repeat_specs_path_)) {
    throw std::invalid_argument("'" + repeat_specs_path_ +
                                "' is not a directory");
  }

  if (!argMap.count("json")) {
    throw std::invalid_argument("--json parameter is required");
  }
  json_path_ = argMap["json"].as<string>();
  DieNotOutputFilePath(json_path_);

  if (!argMap.count("vcf")) {
    throw std::invalid_argument("--vcf parameter is required");
  }
  vcf_path_ = argMap["vcf"].as<string>();
  DieNotOutputFilePath(vcf_path_);

  if (!argMap.count("log")) {
    throw std::invalid_argument("--log parameter is required");
  }
  log_path_ = argMap["log"].as<string>();
  DieNotOutputFilePath(log_path_);

  if (argMap.count("region-extension-length")) {
    region_extension_len_ = argMap["region-extension-length"].as<int>();
  }

  if (argMap.count("min-score")) {
    min_wp_ = argMap["min-score"].as<float>();
    if (min_wp_ > 1) {
      throw std::invalid_argument("min-score must be less than or equal to 1");
    }
  }

  if (argMap.count("min-baseq")) {
    min_baseq_ = argMap["min-baseq"].as<int>();
  }

  if (argMap.count("min-anchor-mapq")) {
    min_anchor_mapq_ = argMap["min-anchor-mapq"].as<int>();
  }

  if (argMap.count("skip-unaligned")) {
    skip_unaligned_ = true;
  }

  if (argMap.count("read-depth")) {
    depth_ = argMap["read-depth"].as<float>();

    if (depth_ < kSmallestPossibleDepth) {
      throw std::invalid_argument("read-depth must be >= " +
                                  lexical_cast<string>(kSmallestPossibleDepth));
    }
  }

  if (argMap.count("sex")) {
    const string sex_encoding = argMap["sex"].as<string>();
    if (sex_encoding == "male") {
      sex_ = Sex::kMale;
    } else if (sex_encoding != "female") {
      throw std::invalid_argument(
          sex_encoding +
          " is an invalid value for sex; it must be either male or female");
    }
  }

  return true;
}
