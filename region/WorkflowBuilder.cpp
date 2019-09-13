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

#include "region/WorkflowBuilder.hh"

#include <memory>
#include <stdexcept>
#include <string>

#include "region/CountingFeature.hh"
#include "region/CountingModel.hh"
#include "region/GraphLocusAnalyzer.hh"
#include "region/StatsAnalyzer.hh"
#include "region/StrAnalyzer.hh"

using std::static_pointer_cast;
using std::string;
using std::unordered_set;
using std::vector;

namespace ehunter
{

LocusAnalyzer::SPtr buildLocusWorkflow(const LocusSpecification& locusSpec, const HeuristicParameters& heuristics)
{
    auto graphLocusPtr = std::make_shared<GraphLocusAnalyzer>();

    if (locusSpec.targetReadExtractionRegions().size() != 1)
    {
        const string message
            = "Locus " + locusSpec.locusId() + " must be associated with exactly one read extraction region";
        throw std::runtime_error(message);
    }

    const auto& regionForGraphModel = locusSpec.targetReadExtractionRegions().front();

    // Construct counting model

    // Initialize counting model
    const int kFlankLength = 1000;
    const int64_t leftFlankEnd = regionForGraphModel.start() + kFlankLength;
    GenomicRegion leftFlank(regionForGraphModel.contigIndex(), regionForGraphModel.start(), leftFlankEnd);

    const int64_t rightFlankStart = regionForGraphModel.end() - kFlankLength;
    GenomicRegion rightFlank(regionForGraphModel.contigIndex(), rightFlankStart, regionForGraphModel.end());

    auto countingModelPtr
        = std::make_shared<CountingModel>(std::initializer_list<GenomicRegion> { leftFlank, rightFlank });

    // TODO: Initialize stats feature
    auto countingFeaturePtr = std::make_shared<CountingFeature>(countingModelPtr);

    auto statsAnalyzer = std::make_shared<StatsAnalyzer>();
    statsAnalyzer->connect(static_pointer_cast<ModelFeature>(countingFeaturePtr));

    // TODO: Add stats analyzer to the locus
    graphLocusPtr->connect(statsAnalyzer);

    // Construct graph model
    auto graphModelPtr = std::make_shared<GraphModel>(regionForGraphModel, locusSpec.regionGraph(), heuristics);

    // Construct graph features and variant analyzers

    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        if (variantSpec.classification().type == VariantType::kRepeat)
        {
            const int motifNode = variantSpec.nodes().front();
            auto strFeaturePtr = std::make_shared<StrFeature>(graphModelPtr, motifNode);
            graphModelPtr->connect(strFeaturePtr.get());

            auto strAnalyzerPtr = std::make_shared<StrAnalyzer>(variantSpec.id());
            strAnalyzerPtr->connect(strFeaturePtr);

            graphLocusPtr->connect(strAnalyzerPtr);

            /*
            weightedPurityCalculators.emplace(std::make_pair(repeatUnit, WeightedPurityCalculator(repeatUnit)));

            if (variantSpec.classification().subtype == VariantSubtype::kRareRepeat)
            {
                if (optionalUnitOfRareRepeat_)
                {
                    const string errorMessage
                        = "Region " + locusSpec_.locusId() + " is not permitted to have more than one rare variant";
                    throw std::logic_error(errorMessage);
                }
                optionalUnitOfRareRepeat_ = repeatUnit;
            } */
        }
        // else if (variantSpec.classification().type == VariantType::kSmallVariant)
        //{
        //    variantAnalyzerPtrs_.emplace_back(new SmallVariantAnalyzer(
        //        variantSpec.id(), variantSpec.classification().subtype, locusSpec.regionGraph(), variantSpec.nodes(),
        //        variantSpec.optionalRefNode(), locusSpec.genotyperParameters()));
        //}
        else
        {
            std::stringstream encoding;
            encoding << variantSpec.classification().type << "/" << variantSpec.classification().subtype;
            throw std::logic_error("Missing logic to create an analyzer for " + encoding.str());
        }
    }

    return graphLocusPtr;
}

vector<RegionModel::SPtr> extractRegionModels(const vector<LocusAnalyzer::SPtr>& locusPtrs)
{
    unordered_set<RegionModel::SPtr> modelPtrs;

    for (const auto& locusPtr : locusPtrs)
    {
        for (const auto& variantPtr : locusPtr->variantAnalyzerPtrs())
        {
            for (const auto& featurePtr : variantPtr->featurePtrs())
            {
                modelPtrs.insert(featurePtr->modelPtr());
            }
        }
    }

    return vector<RegionModel::SPtr>(modelPtrs.begin(), modelPtrs.end());
}

}