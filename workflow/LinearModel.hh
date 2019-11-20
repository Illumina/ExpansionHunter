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

#include "graph_components/ReadClassifier.hh"
#include "workflow/RegionModel.hh"

namespace ehunter
{

class LinearFeature;

class LinearModel : public RegionModel
{
public:
    LinearModel() = delete;
    explicit LinearModel(std::vector<GenomicRegion> readExtractionRegions);
    ~LinearModel() override;

    void addFeature(LinearFeature* feature);
    std::vector<Feature*> modelFeatures() override;

    void analyze(const MappedRead& read, const MappedRead& mate) override;
    void analyze(MappedRead read) override;

private:
    std::vector<LinearFeature*> features_;
    ReadClassifier proximityClassifier_;
};

}
