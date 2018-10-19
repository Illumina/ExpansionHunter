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

#pragma once

#include <fstream>
#include <sstream>
#include <string>

#include "common/genomic_region.h"
#include "region_spec/RegionSpec.hh"

class Outputs
{
public:
    Outputs(const std::string vcf_path, const std::string json_path, const std::string log_path);
    std::ostream& vcf() { return vcf_; }
    std::ostream& json() { return json_; }
    std::ostream& log() { return log_; }

private:
    std::ofstream vcf_;
    std::ofstream json_;
    std::ofstream log_;
};

class Parameters
{
public:
    const double kMinDepthAllowed = 5.0;
    Parameters()
        : region_extension_len_(1000)
        , min_wp_(0.90)
        , min_baseq_(20)
        , min_anchor_mapq_(60)
        , readLen_(-1)
        , haplotype_depth_(0.0)
        , sex_(Sex::kFemale)
        , skip_unaligned_(false)
        , perform_gapped_alignment_(true)
    {
    }
    bool Load(int argc, char** argv);
    std::string bam_path() const { return bam_path_; }
    std::string genome_path() const { return genome_path_; }
    int region_extension_len() const { return region_extension_len_; }
    double min_wp() const { return min_wp_; }
    void set_min_wp(double min_wp) { min_wp_ = min_wp; }
    int min_baseq() const { return min_baseq_; }
    void set_min_baseq(int min_baseq) { min_baseq_ = min_baseq; }
    int min_anchor_mapq() const { return min_anchor_mapq_; }
    bool performGappedAlignment() const { return perform_gapped_alignment_; }
    bool skip_unaligned() const { return skip_unaligned_; }
    int readLen() const { return readLen_; };
    void set_readLen(int readLen) { readLen_ = readLen; }
    double haplotype_depth() const { return haplotype_depth_; }
    void set_haplotype_depth(double haplotype_depth) { haplotype_depth_ = haplotype_depth; }
    std::string sample_name() const { return sample_name_; }
    std::string region_specs_path() const { return region_specs_path_; }
    std::string vcf_path() const { return vcf_path_; }
    std::string json_path() const { return json_path_; }
    std::string log_path() const { return log_path_; }
    bool is_haplotype_depth_set() const { return haplotype_depth_ >= kMinDepthAllowed; }
    Sex sex() const { return sex_; }
    const std::string& alignerName() const { return aligner_name_; }

private:
    std::string bam_path_;
    std::string genome_path_;
    // Maximum distance from a region to search for relevant reads.
    int region_extension_len_;
    double min_wp_;
    int min_baseq_;
    int min_anchor_mapq_;
    int readLen_;
    double haplotype_depth_;
    Sex sex_;
    bool skip_unaligned_;
    bool perform_gapped_alignment_;
    std::string region_specs_path_;
    std::string sample_name_;
    std::string vcf_path_;
    std::string json_path_;
    std::string log_path_;
    std::string aligner_name_;
};
