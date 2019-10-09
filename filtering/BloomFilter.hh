//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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
//

#include <array>
#include <string>
#include <vector>

#include "thirdparty/murmur/MurmurHash3.h"

#include "graphcore/Graph.hh"

namespace ehunter
{

class BloomFilter
{
public:
    BloomFilter();

    using IndexTuple = std::array<uint64_t, 2>;
    IndexTuple computeIndexes(const std::string& kmer) const
    {
        IndexTuple tuple;
        MurmurHash3_x64_128(reinterpret_cast<const unsigned char*>(kmer.data()), kmer.length(), 0, tuple.data());

        tuple[1] = (tuple[0] + tuple[1]) % numBits_;
        tuple[0] = tuple[0] % numBits_;

        return tuple;
    }

    void add(const std::string& kmer)
    {
        auto indexTuple = computeIndexes(kmer);
        bits_[indexTuple[0]] = true;
        bits_[indexTuple[1]] = true;
    }

    bool maybeContains(const std::string& kmer) const
    {
        auto indexTuple = computeIndexes(kmer);
        return bits_[indexTuple[0]] && bits_[indexTuple[1]];
    }

private:
    unsigned int numBits_;
    std::vector<bool> bits_;
};

BloomFilter build(const graphtools::Graph& graph, int kmerLength);

}
