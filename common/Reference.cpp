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

#include "common/Reference.hh"

#include <algorithm>
#include <stdexcept>
#include <string>

using std::string;

namespace ehunter
{

FastaReference::FastaReference(const string& genome_path)
    : genome_path_(genome_path)
{
    fai_ptr_ = fai_load(genome_path_.c_str());
}

FastaReference::~FastaReference() { fai_destroy(fai_ptr_); }

string FastaReference::getSequence(const string& chrom, pos_t start, pos_t end) const
{
    int len; // throwaway...
    // this htslib function is 0-based closed but we need to be open
    char* ref_tmp = faidx_fetch_seq(fai_ptr_, chrom.c_str(), start, end - 1, &len);

    if (!ref_tmp || len < 0 || static_cast<uint>(len) < end - start)
    {
        string region(chrom + ":" + std::to_string(start) + "-" + std::to_string(end));
        throw std::runtime_error(
            "Cannot extract " + region + " from " + genome_path_ + "; chromosome names must match exactly "
            + "(e.g. \"chr1\" and \"1\" are distinct names) "
            + "and coordinates cannot be past the end of the chromosome");
    }
    string sequence("N", len);
    std::transform(ref_tmp, ref_tmp + len, sequence.begin(), ::toupper);
    free(ref_tmp);
    return sequence;
}

}
