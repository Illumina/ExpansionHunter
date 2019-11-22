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

#include "workflow/LinearModel.hh"

#include "spdlog/spdlog.h"
#include "workflow/LinearFeature.hh"

using std::vector;

namespace ehunter
{

LinearModel::LinearModel(std::vector<GenomicRegion> readExtractionRegions)
    : RegionModel(std::move(readExtractionRegions))
    , proximityClassifier_(readExtractionRegions_)
{
}

void LinearModel::analyze(const MappedRead& read, const MappedRead& mate)
{
    analyze(read);
    analyze(mate);
}

void LinearModel::analyze(MappedRead read)
{
    const RegionProximity proximity = proximityClassifier_.classify(read);
    if (proximity == RegionProximity::kInside)
    {
        for (const auto& feature : features_)
        {
            feature->summarize(read);
        }
    }
}

vector<Feature*> LinearModel::modelFeatures()
{
    vector<Feature*> modelFeatures;

    for (const auto& countingFeature : features_)
    {
        modelFeatures.push_back(countingFeature);
    }

    return modelFeatures;
}

LinearModel::~LinearModel()
{
    /*    std::ostringstream regionEncoding;
        regionEncoding << readExtractionRegions_.front();
        regionEncoding << " - ";
        regionEncoding << readExtractionRegions_.back();
        spdlog::info("Region = {}", regionEncoding.str());
        spdlog::info("\tcountReads() = {}", countReads());
        spdlog::info("\tcalculateReadLength() = {}", calculateReadLength());
        spdlog::info("\tcalculateDepth() = {}", calculateDepth()); */
}

void LinearModel::addFeature(LinearFeature* feature) { features_.push_back(feature); }

}
