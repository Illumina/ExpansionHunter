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

#include "region/GraphModel.hh"

#include "alignment/AlignmentFilters.hh"
#include "region/GraphFeature.hh"
#include "spdlog/spdlog.h"

namespace ehunter
{

using graphtools::GappedGraphAligner;
using graphtools::Graph;
using graphtools::GraphAlignment;
using std::list;
using std::string;

GraphModel::GraphModel(GenomicRegion referenceRegion, Graph graph, const HeuristicParameters& heuristics)
    : RegionModel({ referenceRegion }, RegionModel::Type::kTarget)
    , graph_(std::move(graph))
    , aligner_(
          &graph_, heuristics.kmerLenForAlignment(), heuristics.paddingLength(), heuristics.seedAffixTrimLength(),
          heuristics.alignerType())
    , orientationPredictor_(&graph_)
{
}

void GraphModel::analyze(Read read, boost::optional<Read> mate)
{
    ++numPairsProcessed_;
    list<GraphAlignment> readAlignments = align(read);
    list<GraphAlignment> mateAlignments;
    if (mate)
    {
        mateAlignments = align(*mate);
    }

    int numMatchingBases = static_cast<int>(read.sequence().length() / 7.5);
    numMatchingBases = std::max(numMatchingBases, 10);
    LinearAlignmentParameters parameters;
    const int kMinNonRepeatAlignmentScore = numMatchingBases * parameters.matchScore;

    if (!checkIfComesFromGraphLocus(readAlignments, mateAlignments, kMinNonRepeatAlignmentScore))
    {
        return;
    }

    if (!readAlignments.empty() && !mateAlignments.empty())
    {
        for (auto& featurePtr : featurePtrs_)
        {
            featurePtr->process(read, readAlignments, *mate, mateAlignments);
        }
    }
}

list<GraphAlignment> GraphModel::align(Read& read) const
{
    OrientationPrediction predictedOrientation = orientationPredictor_.predict(read.sequence());

    if (predictedOrientation == OrientationPrediction::kAlignsInReverseComplementOrientation)
    {
        read.reverseComplement();
    }
    else if (predictedOrientation == OrientationPrediction::kDoesNotAlign)
    {
        return {};
    }

    return aligner_.align(read.sequence());
}

void GraphModel::addFeature(GraphFeature* featurePtr) { featurePtrs_.push_back(featurePtr); }

std::vector<ModelFeature*> GraphModel::modelFeatures()
{
    std::vector<ModelFeature*> modelFeatures;
    for (const auto& graphFeature : featurePtrs_)
    {
        modelFeatures.push_back(graphFeature);
    }

    return modelFeatures;
}

GraphModel::~GraphModel()
{
    std::ostringstream message;
    message << "Model of " << readExtractionRegions_.front() << " processed " << numPairsProcessed_ << " read pairs";
    spdlog::info(message.str());
}

}