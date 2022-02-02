//
// ExpansionHunter
// Copyright 2016-2021 Illumina, Inc.
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
//

#pragma once

#include <string>
#include <vector>

#include "graphalign/GraphAlignment.hh"

namespace ehunter
{
namespace locus
{

struct AlignmentBufferData
{
    std::string read;
    bool isReversed;
    graphtools::GraphAlignment readAlignment;
};

/// \brief Buffer for read alignments at a single locus
///
/// These read alignments are not needed for standard repeat expansion calling, but are stored for special
/// locus-specific calling extensions.
///
/// This buffer has currently been written to fulfil the requirements of the RFC1 motif analyzer only.
///
class AlignmentBuffer
{
public:
    using buffer_t = std::vector<AlignmentBufferData>;

    /// Test if the given read meets inclusion criteria, and if so, store it in the buffer
    void testAndPushRead(const std::string& read, bool isReversed, const graphtools::GraphAlignment& readAlignment);

    const buffer_t& getBuffer() const { return bufData_; }

private:
    buffer_t bufData_;
};

}
}
