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

#include "reads/Read.hh"

#include <stdexcept>

using std::string;

namespace ehunter
{


namespace reads
{

bool operator==(const Read& read_a, const Read& read_b)
{
    const bool is_id_equal = read_a.read_id == read_b.read_id;
    const bool is_sequence_equal = read_a.sequence == read_b.sequence;
    const bool is_mate_order_equal = read_a.is_first_mate == read_b.is_first_mate;
    return (is_id_equal && is_sequence_equal && is_mate_order_equal);
}

bool operator==(const LinearAlignmentStats& alignment_stats_a, const LinearAlignmentStats& alignment_stats_b)
{
    const bool is_chrom_id_equal = alignment_stats_a.chrom_id == alignment_stats_b.chrom_id;
    const bool is_pos_equal = alignment_stats_a.pos == alignment_stats_b.pos;
    const bool is_mapq_equal = alignment_stats_a.mapq == alignment_stats_b.mapq;
    const bool is_mate_chrom_id_equal = alignment_stats_a.mate_chrom_id == alignment_stats_b.mate_chrom_id;
    const bool is_mate_pos_equal = alignment_stats_a.mate_pos == alignment_stats_b.mate_pos;
    const bool is_mapping_status_equal = alignment_stats_a.is_mapped == alignment_stats_b.is_mapped;
    const bool is_mate_mapping_status_equal = alignment_stats_a.is_mate_mapped == alignment_stats_b.is_mate_mapped;

    return (
        is_chrom_id_equal && is_pos_equal && is_mapq_equal && is_mate_chrom_id_equal && is_mate_pos_equal
        && is_mapping_status_equal && is_mate_mapping_status_equal);
}

std::ostream& operator<<(std::ostream& os, const Read& read)
{
    os << read.read_id << " " << read.sequence;
    return os;
}

} // namespace reads

}
