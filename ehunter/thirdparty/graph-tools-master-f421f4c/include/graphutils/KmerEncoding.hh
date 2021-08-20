//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Chris Saunders <csaunders@illumina.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>

namespace graphtools
{

struct TwoBitKmerEncoder
{
    using KmerKey_t = uint32_t;

    explicit TwoBitKmerEncoder(const size_t kmerLength)
        : kmerLength_(kmerLength)
    {
        const size_t maxKeyBitCount(8 * sizeof(KmerKey_t));
        if (maxKeyBitCount < (kmerLength_ * 2))
        {
            throw std::logic_error(
                "Can't support kmer size of " + std::to_string(kmerLength_) + " with a "
                + std::to_string(maxKeyBitCount) + "bit key type.");
        }
    }

    KmerKey_t encode(const std::string& kmer) const
    {
        if (kmer.size() != kmerLength_)
        {
            throw std::logic_error(
                "kmer size (" + std::to_string(kmer.size()) + ") does not match expected size ("
                + std::to_string(kmerLength_) + "), for kmer '" + kmer + "'.");
        }

        KmerKey_t kmerKey(0);
        for (unsigned i(0); i < kmerLength_; ++i)
        {
            kmerKey = (kmerKey << 2) | baseToIndex(kmer[i]);
        }
        return kmerKey;
    }

    std::string decode(KmerKey_t kmerKey) const
    {
        std::string kmer(kmerLength_, 'N');
        for (unsigned i(0); i < kmerLength_; ++i)
        {
            kmer[kmerLength_ - (i + 1)] = indexToBase(kmerKey & 0x3);
            kmerKey >>= 2;
        }
        return kmer;
    }

private:
    static uint8_t baseToIndex(const char c) { return baseToIndex_.table[static_cast<uint8_t>(c)]; }

    static char indexToBase(const uint8_t i)
    {
        static const char bases[] = "ACGT";
        if (i > 3)
        {
            throw std::logic_error("Unexpected kmer index: '" + std::to_string(i) + "'");
        }
        return bases[i];
    }

    struct BaseToIndex
    {
        BaseToIndex()
        {
            std::fill(table, table + 256, 0);
            table[static_cast<uint8_t>('C')] = 1;
            table[static_cast<uint8_t>('G')] = 2;
            table[static_cast<uint8_t>('T')] = 3;
        }

        uint8_t table[256];
    };

    static BaseToIndex baseToIndex_;
    size_t kmerLength_;
};
}
