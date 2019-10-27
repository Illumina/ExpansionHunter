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

#include "workflow/CnvLocusAnalyzer.hh"
#include "workflow/CnvVariantAnalyzer.hh"
#include "workflow/GraphLocusAnalyzer.hh"
#include "workflow/GraphModel.hh"
#include "workflow/GraphSmallVariant.hh"
#include "workflow/GraphSmallVariantAnalyzer.hh"
#include "workflow/GraphStr.hh"
#include "workflow/GraphStrAnalyzer.hh"
#include "workflow/IrrPairDetector.hh"
#include "workflow/LinearModel.hh"
#include "workflow/ReadCountAnalyzer.hh"
#include "workflow/ReadCounter.hh"

using std::make_shared;
using std::runtime_error;
using std::shared_ptr;
using std::static_pointer_cast;
using std::string;
using std::unordered_set;
using std::vector;

namespace ehunter
{

static shared_ptr<ReadCountAnalyzer>
createStatsAnalyzer(CopyNumberBySex copyNumber, const GenomicRegion& locusLocation, int flankLength)
{
    const int64_t leftFlankStart = locusLocation.start() - flankLength;
    const int64_t leftFlankEnd = locusLocation.start();
    GenomicRegion leftFlank(locusLocation.contigIndex(), leftFlankStart, leftFlankEnd);

    const int64_t rightFlankStart = locusLocation.end();
    const int64_t rightFlankEnd = locusLocation.end() + flankLength;
    GenomicRegion rightFlank(locusLocation.contigIndex(), rightFlankStart, rightFlankEnd);

    vector<GenomicRegion> baselineRegions = { leftFlank, rightFlank };
    auto linearModel = make_shared<LinearModel>(baselineRegions);
    auto readCounter = make_shared<ReadCounter>(linearModel, baselineRegions);
    linearModel->addFeature(readCounter.get());
    return make_shared<ReadCountAnalyzer>(copyNumber, readCounter);
}

static shared_ptr<GraphStrAnalyzer>
createStrAnalyzer(const shared_ptr<GraphModel>& graphModel, const VariantSpecification& variantSpec)
{
    const int motifNode = variantSpec.nodes().front();
    auto str = make_shared<GraphStr>(graphModel, motifNode);
    graphModel->addGraphFeature(str.get());
    auto strAnalyzer = make_shared<GraphStrAnalyzer>(str, variantSpec.id());

    if (variantSpec.classification().subtype == VariantSubtype::kRareRepeat)
    {
        const string& motif = graphModel->graph().nodeSeq(motifNode);
        auto irrPairDetector = make_shared<IrrPairDetector>(graphModel, motif);
        graphModel->addOfftargetReadProcessor(irrPairDetector.get());
        strAnalyzer->addPairedIrrFeature(irrPairDetector);
    }

    return strAnalyzer;
}

static shared_ptr<GraphSmallVariantAnalyzer>
createSmallVariantAnalyzer(const shared_ptr<GraphModel>& graphModel, const VariantSpecification& variantSpec)
{
    auto smallVariant = make_shared<GraphSmallVariant>(graphModel, variantSpec.nodes());
    graphModel->addGraphFeature(smallVariant.get());

    auto smallVariantAnalyzer = make_shared<GraphSmallVariantAnalyzer>(
        smallVariant, variantSpec.id(), variantSpec.classification().subtype, variantSpec.optionalRefNode());

    return smallVariantAnalyzer;
}

shared_ptr<LocusAnalyzer> buildGraphLocusWorkflow(
    const GraphLocusSpecification& locusSpec, const HeuristicParameters& heuristics, BamletWriterPtr bamletWriter)
{
    const auto& locusLocation = locusSpec.locusLocation();

    const double minLocusCoverage = locusSpec.genotyperParameters().minLocusCoverage;
    auto locus = make_shared<GraphLocusAnalyzer>(minLocusCoverage, locusSpec.locusId());
    auto statsAnalyzer
        = createStatsAnalyzer(locusSpec.copyNumberBySex(), locusLocation, heuristics.regionExtensionLength());
    locus->setStats(statsAnalyzer);

    auto graphModel = make_shared<GraphModel>(
        locusSpec.locusId(), locusSpec.targetReadExtractionRegions(), locusSpec.offtargetReadExtractionRegions(),
        locusSpec.regionGraph(), heuristics, bamletWriter);
    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        if (variantSpec.classification().type == VariantType::kRepeat)
        {
            locus->addAnalyzer(createStrAnalyzer(graphModel, variantSpec));
        }
        else if (variantSpec.classification().type == VariantType::kSmallVariant)
        {
            locus->addAnalyzer(createSmallVariantAnalyzer(graphModel, variantSpec));
        }
        else
        {
            std::stringstream encoding;
            encoding << variantSpec.classification().type << "/" << variantSpec.classification().subtype;
            throw runtime_error("Variant " + variantSpec.id() + " is of unknown type " + encoding.str());
        }
    }

    return locus;
}

shared_ptr<LocusAnalyzer> buildCnvLocusWorkflow(
    const CnvLocusSpecification& locusSpec, DepthNormalizer genomeDepthNormalizer,
    const HeuristicParameters& heuristics)
{
    const auto& locusLocation = locusSpec.locusLocation();

    const double minLocusCoverage = locusSpec.genotyperParameters().minLocusCoverage;
    auto locus = make_shared<CnvLocusAnalyzer>(minLocusCoverage, locusSpec.locusId(), locusSpec.locusSubtype());
    auto statsAnalyzer
        = createStatsAnalyzer(locusSpec.copyNumberBySex(), locusLocation, heuristics.regionExtensionLength());
    locus->setStats(statsAnalyzer);

    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        if (variantSpec.classification().type == VariantType::kCNV)
        {
            double regionLength = variantSpec.referenceLocus().end() - variantSpec.referenceLocus().start();

            assert(variantSpec.parameters());
            CnvGenotyperParameters cnvParameters = *variantSpec.parameters();

            vector<GenomicRegion> variantRegion { variantSpec.referenceLocus() };
            auto linearModel = make_shared<LinearModel>(variantRegion);
            auto readCounter = make_shared<ReadCounter>(linearModel, variantRegion);
            linearModel->addFeature(readCounter.get());
            locus->addAnalyzer(make_shared<CnvVariantAnalyzer>(
                variantSpec.id(), regionLength, variantSpec.classification().subtype, locusSpec.copyNumberBySex(),
                cnvParameters, readCounter, genomeDepthNormalizer));
        }
        else
        {
            std::stringstream encoding;
            encoding << variantSpec.classification().type << "/" << variantSpec.classification().subtype;
            throw runtime_error("Variant " + variantSpec.id() + " is of unknown type " + encoding.str());
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
