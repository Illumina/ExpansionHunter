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

#include "workflow/GraphLocusAnalyzer.hh"
#include "workflow/GraphModel.hh"
#include "workflow/GraphSmallVariant.hh"
#include "workflow/GraphSmallVariantAnalyzer.hh"
#include "workflow/GraphStr.hh"
#include "workflow/GraphStrAnalyzer.hh"
#include "workflow/LinearModel.hh"
#include "workflow/LinearModelFeature.hh"
#include "workflow/OfftargetFeature.hh"
#include "workflow/StatsAnalyzer.hh"

using std::shared_ptr;
using std::static_pointer_cast;
using std::string;
using std::unordered_set;
using std::vector;

namespace ehunter
{

shared_ptr<LocusAnalyzer> buildLocusWorkflow(const LocusSpecification& locusSpec, const HeuristicParameters& heuristics)
{
    if (locusSpec.targetReadExtractionRegions().size() != 1)
    {
        const string message = "Locus " + locusSpec.locusId() + " must be associated with one read extraction region";
        throw std::runtime_error(message);
    }

    const auto& extractionRegion = locusSpec.targetReadExtractionRegions().front();

    // Initialize counting model
    const int64_t leftFlankEnd = extractionRegion.start() + heuristics.regionExtensionLength();
    GenomicRegion leftFlank(extractionRegion.contigIndex(), extractionRegion.start(), leftFlankEnd);

    const int64_t rightFlankStart = extractionRegion.end() - heuristics.regionExtensionLength();
    GenomicRegion rightFlank(extractionRegion.contigIndex(), rightFlankStart, extractionRegion.end());

    // Initialize stats feature
    vector<GenomicRegion> baselineRegions = { leftFlank, rightFlank };
    auto linearModel = std::make_shared<LinearModel>(baselineRegions);
    auto countingFeature = std::make_shared<LinearModelFeature>(linearModel, baselineRegions);
    linearModel->addFeature(countingFeature.get());
    auto statsAnalyzer = std::make_shared<StatsAnalyzer>(countingFeature);

    auto locus = std::make_shared<GraphLocusAnalyzer>(locusSpec.locusId());
    locus->setStats(statsAnalyzer);

    auto graphModel = std::make_shared<GraphModel>(extractionRegion, locusSpec.regionGraph(), heuristics);

    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        if (variantSpec.classification().type == VariantType::kRepeat)
        {
            const int motifNode = variantSpec.nodes().front();
            auto strFeaturePtr = std::make_shared<GraphStr>(graphModel, motifNode);
            graphModel->addVariant(strFeaturePtr.get());

            auto strAnalyzerPtr = std::make_shared<GraphStrAnalyzer>(strFeaturePtr, variantSpec.id());
            locus->addAnalyzer(strAnalyzerPtr);

            if (variantSpec.classification().subtype == VariantSubtype::kRareRepeat)
            {
                const string& motif = graphModel->graph().nodeSeq(motifNode);
                auto pairedIrrFeaturePtr = std::make_shared<OfftargetFeature>(graphModel, motif);
                graphModel->addOfftargetFeature(pairedIrrFeaturePtr.get());
                strAnalyzerPtr->addPairedIrrFeature(pairedIrrFeaturePtr);
            }
        }
        else if (variantSpec.classification().type == VariantType::kSmallVariant)
        {
            auto smallVariantFeature = std::make_shared<GraphSmallVariant>(graphModel, variantSpec.nodes());
            graphModel->addVariant(smallVariantFeature.get());

            auto smallVariantAnalyzer = std::make_shared<GraphSmallVariantAnalyzer>(
                smallVariantFeature, variantSpec.id(), variantSpec.classification().subtype,
                variantSpec.optionalRefNode());
        }
        else
        {
            std::stringstream encoding;
            encoding << variantSpec.classification().type << "/" << variantSpec.classification().subtype;
            throw std::logic_error("Missing logic to create an analyzer for " + encoding.str());
        }
    }

    return locus;
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