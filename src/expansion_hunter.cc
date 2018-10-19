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

#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/filesystem.hpp>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "thirdparty/spdlog/fmt/ostr.h"
#include "thirdparty/spdlog/spdlog.h"

#include "common/parameters.h"
#include "output/JsonWriter.hh"
#include "output/VcfWriter.hh"
#include "region_analysis/RepeatFindings.hh"
#include "sample_analysis/HtsSeekingSampleAnalyzer.hh"
#include "sample_analysis/HtsStreamingSampleAnalyzer.hh"
#include "src/version.h"

namespace spd = spdlog;

using std::map;
using std::string;
using std::vector;

// Returns the length of the first read in a BAM file.
size_t extractReadLength(const string& bam_path)
{
    // Open a BAM file for reading.
    samFile* file_ptr = sam_open(bam_path.c_str(), "r");
    if (!file_ptr)
    {
        throw std::runtime_error("Failed to read BAM file '" + bam_path + "'");
    }
    bam_hdr_t* header_ptr = sam_hdr_read(file_ptr);
    if (!header_ptr)
    {
        throw std::runtime_error("BamFile::Init: Failed to read BAM header: '" + bam_path + "'");
    }

    enum
    {
        kSupplimentaryAlign = 0x800,
        kSecondaryAlign = 0x100
    };

    size_t readLen = 99;
    bam1_t* align_ptr = bam_init1();
    int ret;
    while ((ret = sam_read1(file_ptr, header_ptr, align_ptr)) >= 0)
    {
        const bool is_supplimentary = align_ptr->core.flag & kSupplimentaryAlign;
        const bool is_secondary = align_ptr->core.flag & kSecondaryAlign;
        const bool is_primary_align = (!is_supplimentary) && (!is_secondary);
        if (is_primary_align)
        {
            readLen = align_ptr->core.l_qseq;
            break;
        }
    }

    if (ret < 0)
    {
        throw std::runtime_error("Failed to extract a read from BAM file");
    }

    bam_destroy1(align_ptr);
    bam_hdr_destroy(header_ptr);
    sam_close(file_ptr);

    return readLen;
}

bool isBamFile(const string& htsFilePath)
{
    string extension = boost::filesystem::extension(htsFilePath);
    return extension == ".bam";
}

int main(int argc, char** argv)
{
    auto console = spd::stderr_color_mt("console");
    spd::set_pattern("%Y-%m-%dT%H:%M:%S,[%v]");

    try
    {
        Parameters parameters;
        console->info("Starting {}", kProgramVersion);

        if (!parameters.Load(argc, argv))
        {
            return 1;
        }

        console->info("Analyzing sample {}", parameters.sample_name());

        Outputs outputs(parameters.vcf_path(), parameters.json_path(), parameters.log_path());

        const int32_t readLength = extractReadLength(parameters.bam_path());
        parameters.set_readLen(readLength);
        console->info("Read length is set to {}", parameters.readLen());

        console->info("Loading repeat specifications from disk {}", parameters.region_specs_path());
        RefGenome reference(parameters.genome_path());
        const RegionCatalog regionSpecs
            = loadRegionSpecsFromDisk(parameters.region_specs_path(), reference, parameters.sex());

        console->info("Running sample analysis");
        SampleFindings sampleFindings;

        if (isBamFile(parameters.bam_path()))
        {
            sampleFindings = htsSeekingSampleAnalysis(
                parameters.bam_path(), parameters.haplotype_depth(), parameters.readLen(), regionSpecs,
                parameters.alignerName(), outputs.log());
        }
        else
        {
            sampleFindings = htslibStreamingSampleAnalyzer(
                parameters.bam_path(), parameters.haplotype_depth(), parameters.readLen(), regionSpecs,
                parameters.alignerName(), outputs.log());
        }

        console->info("Writing output to disk");
        VcfWriter vcfWriter(parameters.sample_name(), parameters.readLen(), regionSpecs, sampleFindings);
        outputs.vcf() << vcfWriter;

        JsonWriter jsonWriter(parameters.sample_name(), parameters.readLen(), regionSpecs, sampleFindings);
        outputs.json() << jsonWriter;
    }
    catch (const std::exception& e)
    {
        console->error(e.what());
        return 1;
    }

    return 0;
}
