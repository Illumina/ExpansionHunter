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
createStatsAnalyzer(CopyNumberBySex copyNumber, const vector<GenomicRegion>& statsRegions)
{
    auto linearModel = make_shared<LinearModel>(statsRegions);
    auto readCounter = make_shared<ReadCounter>(linearModel, statsRegions);
    linearModel->addFeature(readCounter.get());
    return make_shared<ReadCountAnalyzer>(copyNumber, readCounter);
}

static shared_ptr<GraphStrAnalyzer>
createStrAnalyzer(const shared_ptr<GraphModel>& graphModel, const GraphVariantSpec& variantSpec)
{
    const int motifNode = variantSpec.nodes().front();
    auto str = make_shared<GraphStr>(graphModel, motifNode);
    graphModel->addGraphFeature(str.get());
    auto strAnalyzer = make_shared<GraphStrAnalyzer>(str, variantSpec.id());

    if (variantSpec.classification().subtype == GraphVariantClassification::Subtype::kRareRepeat)
    {
        const string& motif = graphModel->graph().nodeSeq(motifNode);
        auto irrPairDetector = make_shared<IrrPairDetector>(graphModel, motif);
        graphModel->addOfftargetReadProcessor(irrPairDetector.get());
        strAnalyzer->addPairedIrrFeature(irrPairDetector);
    }

    return strAnalyzer;
}

static shared_ptr<GraphSmallVariantAnalyzer>
createSmallVariantAnalyzer(const shared_ptr<GraphModel>& graphModel, const GraphVariantSpec& variantSpec)
{
    auto smallVariant = make_shared<GraphSmallVariant>(graphModel, variantSpec.nodes());
    graphModel->addGraphFeature(smallVariant.get());

    auto smallVariantAnalyzer = make_shared<GraphSmallVariantAnalyzer>(
        smallVariant, variantSpec.id(), variantSpec.classification().subtype, variantSpec.optionalRefNode());

    return smallVariantAnalyzer;
}

shared_ptr<LocusAnalyzer> buildGraphLocusWorkflow(
    const GraphLocusSpec& locusSpec, const HeuristicParameters& heuristics, BamletWriterPtr bamletWriter)
{
    const double minLocusCoverage = locusSpec.genotyperParams().minLocusCoverage;
    auto locus = make_shared<GraphLocusAnalyzer>(minLocusCoverage, locusSpec.locusId());
    auto statsAnalyzer = createStatsAnalyzer(locusSpec.copyNumberBySex(), locusSpec.analysisRegions().statsRegions);
    locus->setStats(statsAnalyzer);

    auto graphModel = make_shared<GraphModel>(
        locusSpec.locusId(), locusSpec.analysisRegions().regionsWithReads,
        locusSpec.analysisRegions().offtargetRegionsWithReads, locusSpec.graph(), heuristics, bamletWriter);
    for (const auto& variantSpec : locusSpec.variants())
    {
        if (variantSpec.classification().type == GraphVariantClassification::Type::kRepeat)
        {
            locus->addAnalyzer(createStrAnalyzer(graphModel, variantSpec));
        }
        else if (variantSpec.classification().type == GraphVariantClassification::Type::kSmallVariant)
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

shared_ptr<LocusAnalyzer> buildCnvLocusWorkflow(const CnvLocusSpec& locusSpec)
{
    auto locus = make_shared<CnvLocusAnalyzer>(locusSpec.locusId(), locusSpec.locusType(), locusSpec.outputVariant());
    auto statsAnalyzer = createStatsAnalyzer(locusSpec.copyNumberBySex(), locusSpec.regionsWithReads());
    locus->setStats(statsAnalyzer);

    for (const auto& variantSpec : locusSpec.variants())
    {

        double regionLength = variantSpec.location().end() - variantSpec.location().start();

        CnvGenotyperParameters cnvParameters = variantSpec.genotyperParams();

        vector<GenomicRegion> variantRegion { variantSpec.location() };
        auto linearModel = make_shared<LinearModel>(variantRegion);
        auto readCounter = make_shared<ReadCounter>(linearModel, variantRegion);
        linearModel->addFeature(readCounter.get());
        locus->addAnalyzer(make_shared<CnvVariantAnalyzer>(
            variantSpec.id(), regionLength, variantSpec.variantType(), locusSpec.copyNumberBySex(),
            cnvParameters, readCounter));
        
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
