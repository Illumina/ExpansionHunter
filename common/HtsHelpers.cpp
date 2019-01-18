//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "common/HtsHelpers.hh"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>

#include "thirdparty/spdlog/spdlog.h"

using std::pair;
using std::string;
using std::vector;

namespace spd = spdlog;

namespace ehunter
{
namespace htshelpers
{

    string decodeQuals(bam1_t* htsAlignPtr)
    {
        string quals;
        uint8_t* htsQualsPtr = bam_get_qual(htsAlignPtr);
        const int readLength = htsAlignPtr->core.l_qseq;
        quals.resize(readLength);

        for (int index = 0; index < readLength; ++index)
        {
            quals[index] = static_cast<char>(33 + htsQualsPtr[index]);
        }

        return quals;
    }

    string decodeBases(bam1_t* htsAlignPtr)
    {
        string bases;
        uint8_t* htsSeqPtr = bam_get_seq(htsAlignPtr);
        const int32_t readLength = htsAlignPtr->core.l_qseq;
        bases.resize(readLength);

        for (int32_t index = 0; index < readLength; ++index)
        {
            bases[index] = seq_nt16_str[bam_seqi(htsSeqPtr, index)];
        }

        return bases;
    }

    LinearAlignmentStats decodeAlignmentStats(bam1_t* htsAlignPtr)
    {
        LinearAlignmentStats alignmentStats;
        alignmentStats.chromId = htsAlignPtr->core.tid;
        alignmentStats.pos = htsAlignPtr->core.pos;
        alignmentStats.mapq = htsAlignPtr->core.qual;
        alignmentStats.mateChromId = htsAlignPtr->core.mtid;
        alignmentStats.matePos = htsAlignPtr->core.mpos;

        uint32_t samFlag = htsAlignPtr->core.flag;
        alignmentStats.isMapped = !(samFlag & SamFlags::kIsUnmapped);
        alignmentStats.isMateMapped = !(samFlag & SamFlags::kIsMateUnmapped);

        return alignmentStats;
    }

    static string lowercaseLowQualityBases(const string& bases, const string& quals, int lowBaseQualityCutoff = 20)
    {
        const int qualityScoreOffset = 33;
        string query = bases;
        for (int index = 0; index != static_cast<int>(bases.size()); ++index)
        {
            if (quals[index] - qualityScoreOffset <= lowBaseQualityCutoff)
            {
                query[index] = std::tolower(bases[index]);
            }
        }
        return query;
    }

    Read decodeRead(bam1_t* htsAlignPtr)
    {
        const uint32_t samFlag = htsAlignPtr->core.flag;
        const bool isFirstMate = samFlag & SamFlags::kIsFirstMate;

        const string fragmentId = bam_get_qname(htsAlignPtr);
        MateNumber mateNumber = isFirstMate ? MateNumber::kFirstMate : MateNumber::kSecondMate;
        ReadId readId(fragmentId, mateNumber);

        string bases = decodeBases(htsAlignPtr);
        string quals = decodeQuals(htsAlignPtr);
        string sequence = lowercaseLowQualityBases(bases, quals);

        return Read(readId, sequence);
    }

    ReferenceContigInfo decodeContigInfo(bam_hdr_t* htsHeaderPtr)
    {
        vector<pair<string, int64_t>> contigNamesAndSizes;
        const int32_t numContigs = htsHeaderPtr->n_targets;
        contigNamesAndSizes.reserve(numContigs);

        for (int32_t contigIndex = 0; contigIndex != numContigs; ++contigIndex)
        {
            const string contig = htsHeaderPtr->target_name[contigIndex];
            int64_t size = htsHeaderPtr->target_len[contigIndex];
            contigNamesAndSizes.push_back(std::make_pair(contig, size));
        }

        return ReferenceContigInfo(contigNamesAndSizes);
    }

} // namespace htshelpers
}
