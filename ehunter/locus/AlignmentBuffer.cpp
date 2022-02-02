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

#include "locus/AlignmentBuffer.hh"

namespace ehunter
{
namespace locus
{

/// \brief Return true if alignment overlaps with segment 1 in the LocusStructure
///
/// Note this isn't generalized to recognize all segments which are repeats in more complex locus structures, it should
/// only work with a single expansion like RFC1.
///
static bool testRepeatMotifOverlap(const graphtools::GraphAlignment& readAlign)
{
    for (size_t index(0); index < readAlign.size(); ++index)
    {
        if (readAlign.path().getNodeIdByIndex(index) == 1)
        {
            return true;
        }
    }
    return false;
}

void AlignmentBuffer::testAndPushRead(
    const std::string& read, const bool isReversed, const graphtools::GraphAlignment& readAlignment)
{
    // Test if the read should be included in the buffer, for now just check that the read touches the repeat at all
    if (not testRepeatMotifOverlap(readAlignment))
    {
        return;
    }

    bufData_.emplace_back(AlignmentBufferData { read, isReversed, readAlignment });
}

}
}
