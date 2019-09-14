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

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/optional.hpp>

#include "region/RegionModel.hh"

namespace ehunter
{

class CountingFeature;

class CountingModel : public RegionModel
{
public:
    CountingModel() = delete;
    explicit CountingModel(std::vector<GenomicRegion> readExtractionRegions);
    ~CountingModel() override = default;

    std::vector<ModelFeature*> modelFeatures() override;

    void analyze(Read read, boost::optional<Read> mate) override;
    int readCount() const;
    int meanReadLength() const;
    double depth() const;

private:
    std::vector<CountingFeature*> featurePtrs_;

    using AccumulatorStats
        = boost::accumulators::features<boost::accumulators::tag::count, boost::accumulators::tag::mean>;
    using Accumulator = boost::accumulators::accumulator_set<int, AccumulatorStats>;

    Accumulator readLengthAccumulator_;
};

}
