//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>
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

#include <memory>
#include <string>
#include <vector>

#include "graphalign/GraphAlignment.hh"
#include "graphcore/GraphReferenceMapping.hh"

namespace graphtools
{

class AlignmentWriter
{
public:
    virtual ~AlignmentWriter() = default;
    virtual void write(
        const std::string& locusId, const std::string& fragmentName, const std::string& query, bool isFirstMate,
        bool isReversed, bool isMateReversed, const GraphAlignment& alignment)
        = 0;
};

class BlankAlignmentWriter : public AlignmentWriter
{
public:
    ~BlankAlignmentWriter() override = default;
    void write(
        const std::string& locusId, const std::string& fragmentName, const std::string& query, bool isFirstMate,
        bool isReversed, bool isMateReversed, const GraphAlignment& alignment) override;
};
}
