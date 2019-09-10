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

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "common/GenomicRegion.hh"
#include "reads/Read.hh"

namespace ehunter
{

class Region
{
public:
    using SPtr = std::shared_ptr<Region>;
    enum class Type
    {
        kTarget,
        kOfftarget
    };

    // Add read extraction regions
    explicit Region(std::string regionId, Type type);
    virtual ~Region() = default;

    Type type() const { return type_; }
    const std::vector<GenomicRegion>& readExtractionRegions() const { return readExtractionRegions_; }

    virtual void analyze(Read read, boost::optional<Read> mate) = 0;

protected:
    std::string regionId_;
    Type type_;
    std::vector<GenomicRegion> readExtractionRegions_;
    // ChromType typeOfChromLocusLocatedOn_;
};

}