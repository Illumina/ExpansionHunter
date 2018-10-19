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

#include <iomanip>
#include <ios>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

using std::string;
using std::vector;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

Outputs::Outputs(const string vcf_path, const string json_path, const string log_path)
{
    vcf_.open(vcf_path.c_str());
    if (!vcf_.is_open())
    {
        throw std::runtime_error("ERROR: Failed to open " + vcf_path + " for writing: " + strerror(errno));
    }

    json_.open(json_path.c_str());
    if (!json_.is_open())
    {
        throw std::runtime_error("ERROR: Failed to open " + json_path + " for writing: " + strerror(errno));
    }

    log_.open(log_path.c_str());
    if (!log_.is_open())
    {
        throw std::runtime_error("ERROR: Failed to open " + log_path + " for writing: " + strerror(errno));
    }
}

static bool CheckIfIndexFileExists(const string& bam_path)
{
    vector<string> kPossibleIndexExtensions = { ".bai", ".csi", ".crai" };

    for (const string& index_extension : kPossibleIndexExtensions)
    {
        if (fs::exists(bam_path + index_extension))
        {
            return true;
        }
    }

    return false;
}

void DieIfOutputPathDoesntExist(const string& output_path_str)
{
    const fs::path output_path(output_path_str);
    const fs::path output_dir(output_path.parent_path());

    const bool is_no_dir = output_dir.empty();
    const bool is_existing_dir = fs::is_directory(output_dir);
    const bool is_valid_fname = fs::portable_posix_name(output_path.filename().string());

    if ((is_no_dir || is_existing_dir) && is_valid_fname)
    {
        return;
    }

    throw std::invalid_argument("ERROR: " + output_path_str + " is not a valid output path");
}

bool Parameters::Load(int argc, char** argv)
{
    // clang-format off
  po::options_description usage("Allowed options");
  usage.add_options()
      ("help", "Print help message")
      ("version", "Print version number")
      ("bam", po::value<string>()->required(), "BAM file")
      ("ref-fasta", po::value<string>()->required(), "FASTA file with reference genome")
      ("repeat-specs", po::value<string>()->required(), "JSON file with repeats to genotype")
      ("vcf", po::value<string>()->required(), "Output VCF file")
      ("json", po::value<string>()->required(), "Output JSON file")
      ("log", po::value<string>()->required(), "Output read alignment file")
      ("region-extension-length", po::value<int>()->default_value(1000), "How far from on/off-target regions to search for informative reads")
      ("min-score", po::value<double>()->default_value(0.90, "0.90"), "Minimum weighted purity score required to flag a read as an in-repeat read; must be between 0 and 1")
      ("min-baseq", po::value<int>()->default_value(20), "Minimum quality of a high-confidence base call")
      ("min-anchor-mapq", po::value<int>()->default_value(60), "Minimum MAPQ of a read anchor")
      ("skip-unaligned", po::bool_switch()->default_value(false), "Skip unaligned reads when searching for IRRs")
      ("read-depth", po::value<double>()->required(), "Read depth")
      ("sex", po::value<string>()->default_value("female"), "Sex of the sample; must be either male or female")
      ("aligner", po::value<string>(&aligner_name_)->default_value("dag-aligner"), "dag-aligner or path-aligner");
    // clang-format on

    if (argc == 1)
    {
        std::cerr << usage << std::endl;
        return false;
    }

    po::variables_map arg_map;
    po::store(po::command_line_parser(argc, argv).options(usage).run(), arg_map);

    if (arg_map.count("help"))
    {
        std::cerr << usage << std::endl;
        return false;
    }

    if (arg_map.count("version"))
    {
        return false;
    }

    po::notify(arg_map);

    bam_path_ = arg_map["bam"].as<string>();
    if (!fs::exists(bam_path_))
    {
        throw std::invalid_argument("ERROR: " + bam_path_ + " does not exist");
    }
    if (!CheckIfIndexFileExists(bam_path_))
    {
        throw std::invalid_argument("ERROR: Could not find index file for BAM: " + bam_path_);
    }

    // Extract sample name.
    fs::path boost_bam_path(bam_path_);
    sample_name_ = boost_bam_path.stem().string();

    genome_path_ = arg_map["ref-fasta"].as<string>();
    if (!fs::exists(genome_path_))
    {
        throw std::invalid_argument("ERROR: " + genome_path_ + " does not exist");
    }

    region_specs_path_ = arg_map["repeat-specs"].as<string>();
    DieIfOutputPathDoesntExist(region_specs_path_);

    json_path_ = arg_map["json"].as<string>();
    DieIfOutputPathDoesntExist(json_path_);

    vcf_path_ = arg_map["vcf"].as<string>();
    DieIfOutputPathDoesntExist(vcf_path_);

    log_path_ = arg_map["log"].as<string>();
    DieIfOutputPathDoesntExist(log_path_);

    region_extension_len_ = arg_map["region-extension-length"].as<int>();

    min_wp_ = arg_map["min-score"].as<double>();
    if (min_wp_ > 1)
    {
        throw std::invalid_argument("min-score must be less than or equal to 1");
    }

    min_baseq_ = arg_map["min-baseq"].as<int>();
    min_anchor_mapq_ = arg_map["min-anchor-mapq"].as<int>();
    skip_unaligned_ = arg_map["skip-unaligned"].as<bool>();

    haplotype_depth_ = arg_map["read-depth"].as<double>() / 2;
    if (haplotype_depth_ * 2 < kMinDepthAllowed)
    {
        throw std::invalid_argument("read-depth must be at least " + std::to_string(kMinDepthAllowed));
    }

    const string sex_encoding = arg_map["sex"].as<string>();
    if (sex_encoding == "male")
    {
        sex_ = Sex::kMale;
    }
    else if (sex_encoding == "female")
    {
        sex_ = Sex::kFemale;
    }
    else
    {
        throw std::invalid_argument("ERROR: " + sex_encoding + " is invalid for sex; must be either male or female");
    }

    return true;
}
