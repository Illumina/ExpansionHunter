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

#include "region/GraphLocus.hh"
#include "region/StrAnalyzer.hh"

using std::string;
using std::unordered_set;
using std::vector;

namespace ehunter
{

LocusAnalyzer::SPtr buildLocusWorkflow(const LocusSpecification& locusSpec, const HeuristicParameters& heuristics)
{
    // TODO: Build GraphLocus
    // TODO: It should be initialized with LocusStats and analyzer ptrs
    auto graphLocus = std::make_shared<GraphLocus>();

    // Construct graph model
    if (locusSpec.targetReadExtractionRegions().size() != 1)
    {
        const string message
            = "Locus " + locusSpec.locusId() + " must be associated with exactly one read extraction region";
        throw std::runtime_error(message);
    }

    const auto& regionForGraphModel = locusSpec.targetReadExtractionRegions().front();
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

            graphLocus->connect(strAnalyzerPtr);

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

    return graphLocus;
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