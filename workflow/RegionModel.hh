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

#include <vector>

#include <boost/optional.hpp>

#include "common/GenomicRegion.hh"
#include "reads/Read.hh"

namespace ehunter
{

class Feature;

class RegionModel
{
public:
    explicit RegionModel(std::vector<GenomicRegion> readExtractionRegions);
    virtual ~RegionModel() = default;

    const std::vector<GenomicRegion>& readExtractionRegions() const { return readExtractionRegions_; }

    virtual void analyze(MappedRead read, MappedRead mate) = 0;
    virtual void analyze(MappedRead read) = 0;

    virtual std::vector<Feature*> modelFeatures() = 0;

protected:
    std::vector<GenomicRegion> readExtractionRegions_;

    // ChromType typeOfChromLocusLocatedOn_;
};

}
