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

#include "workflow/CountingModel.hh"

#include "spdlog/spdlog.h"
#include "workflow/CountingFeature.hh"

using std::vector;

namespace ehunter
{

CountingModel::CountingModel(std::vector<GenomicRegion> readExtractionRegions)
    : RegionModel(std::move(readExtractionRegions), Type::kTarget)
    , proximityClassifier_(readExtractionRegions_)
{
}

void CountingModel::analyze(MappedRead read, MappedRead mate)
{
    analyze(read);
    analyze(mate);
}

void CountingModel::analyze(MappedRead read)
{
    const RegionProximity proximity = proximityClassifier_.classify(read);
    if (proximity == RegionProximity::kInside)
    {
        for (const auto& feature : featurePtrs_)
        {
            feature->addReadInfo(read.sequence().length());
        }
    }
}

vector<ModelFeature*> CountingModel::modelFeatures()
{
    vector<ModelFeature*> modelFeatures;

    for (const auto& countingFeature : featurePtrs_)
    {
        modelFeatures.push_back(countingFeature);
    }

    return modelFeatures;
}

CountingModel::~CountingModel()
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

void CountingModel::addFeature(CountingFeature* featurePtr) { featurePtrs_.push_back(featurePtr); }

}
