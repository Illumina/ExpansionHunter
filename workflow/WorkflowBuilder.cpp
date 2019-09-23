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

#include "workflow/WorkflowBuilder.hh"

#include <stdexcept>
#include <string>

#include "workflow/CountingFeature.hh"
#include "workflow/CountingModel.hh"
#include "workflow/GraphLocusAnalyzer.hh"
#include "workflow/GraphModel.hh"
#include "workflow/PairedIrrFeature.hh"
#include "workflow/SmallVariantAnalyzer.hh"
#include "workflow/SmallVariantFeature.hh"
#include "workflow/StatsAnalyzer.hh"
#include "workflow/StrAnalyzer.hh"
#include "workflow/StrFeature.hh"

using std::shared_ptr;
using std::static_pointer_cast;
using std::string;
using std::unordered_set;
using std::vector;

namespace ehunter
{

shared_ptr<LocusAnalyzer> buildLocusWorkflow(const LocusSpecification& locusSpec, const HeuristicParameters& heuristics)
{
    auto graphLocusPtr = std::make_shared<GraphLocusAnalyzer>(locusSpec.locusId());

    if (locusSpec.targetReadExtractionRegions().size() != 1)
    {
        const string message
            = "Locus " + locusSpec.locusId() + " must be associated with exactly one read extraction workflow";
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

    // Initialize stats feature
    vector<GenomicRegion> statsRegions = { leftFlank, rightFlank };
    auto countingModelPtr = std::make_shared<CountingModel>(statsRegions);
    auto countingFeaturePtr = std::make_shared<CountingFeature>(countingModelPtr, statsRegions);
    countingModelPtr->addFeature(countingFeaturePtr.get());
    auto statsAnalyzerPtr = std::make_shared<StatsAnalyzer>(countingFeaturePtr);
    graphLocusPtr->setStats(statsAnalyzerPtr);

    // Construct graph model
    auto graphModelPtr = std::make_shared<GraphModel>(regionForGraphModel, locusSpec.regionGraph(), heuristics);

    // Construct graph features and variant analyzers

    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        if (variantSpec.classification().type == VariantType::kRepeat)
        {
            const int motifNode = variantSpec.nodes().front();
            auto strFeaturePtr = std::make_shared<StrFeature>(graphModelPtr, motifNode);
            graphModelPtr->addFeature(strFeaturePtr.get());

            auto strAnalyzerPtr = std::make_shared<StrAnalyzer>(strFeaturePtr, variantSpec.id());
            graphLocusPtr->addAnalyzer(strAnalyzerPtr);

            if (variantSpec.classification().subtype == VariantSubtype::kRareRepeat)
            {
                const string& motif = graphModelPtr->graph().nodeSeq(motifNode);
                auto pairedIrrFeaturePtr = std::make_shared<PairedIrrFeature>(graphModelPtr, motif);
                graphModelPtr->addPairedIrrFeature(pairedIrrFeaturePtr.get());
                strAnalyzerPtr->addPairedIrrFeature(pairedIrrFeaturePtr);
            }
        }
        else if (variantSpec.classification().type == VariantType::kSmallVariant)
        {
            auto smallVariantFeature = std::make_shared<SmallVariantFeature>(graphModelPtr, variantSpec.nodes());
            graphModelPtr->addFeature(smallVariantFeature.get());

            auto smallVariantAnalyzer = std::make_shared<SmallVariantAnalyzer>(
                smallVariantFeature, variantSpec.id(), variantSpec.classification().subtype,
                variantSpec.optionalRefNode());
            //    variantAnalyzerPtrs_.emplace_back(new SmallVariantAnalyzer(
            //        variantSpec.id(), variantSpec.classification().subtype, locusSpec.regionGraph(),
            //        variantSpec.nodes(), variantSpec.optionalRefNode(), locusSpec.genotyperParameters()));
        }
        else
        {
            std::stringstream encoding;
            encoding << variantSpec.classification().type << "/" << variantSpec.classification().subtype;
            throw std::logic_error("Missing logic to create an analyzer for " + encoding.str());
        }
    }

    return graphLocusPtr;
}

vector<shared_ptr<RegionModel>> extractRegionModels(const vector<shared_ptr<LocusAnalyzer>>& loci)
{
    unordered_set<shared_ptr<RegionModel>> models;

    for (const auto& locus : loci)
    {
        for (const auto& featureAnalyzer : locus->featureAnalyzers())
        {
            for (const auto& feature : featureAnalyzer->features())
            {
                models.insert(feature->model());
            }
        }
    }

    return vector<shared_ptr<RegionModel>>(models.begin(), models.end());
}
}