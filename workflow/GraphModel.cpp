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

#include "workflow/GraphModel.hh"

#include "alignment/AlignmentFilters.hh"
#include "spdlog/spdlog.h"
#include "workflow/GraphFeature.hh"
#include "workflow/IrrPairDetector.hh"

namespace ehunter
{

using graphtools::GappedGraphAligner;
using graphtools::Graph;
using graphtools::GraphAlignment;
using std::list;
using std::string;

GraphModel::GraphModel(GenomicRegion referenceRegion, Graph graph, const HeuristicParameters& heuristics)
    : RegionModel({ referenceRegion })
    , readClassifier_(readExtractionRegions_)
    , graph_(std::move(graph))
    , aligner_(
          &graph_, heuristics.kmerLenForAlignment(), heuristics.paddingLength(), heuristics.seedAffixTrimLength(),
          heuristics.alignerType())
    , orientationPredictor_(&graph_)
{
}

void GraphModel::analyze(MappedRead read, MappedRead mate)
{
    RegionProximity type = readClassifier_.classify(read, mate);

    if (type == RegionProximity::kOverlapsOrNear)
    {
        return;
    }

    if (type == RegionProximity::kFar && offtargetProcessor_ != nullptr)
    {
        offtargetProcessor_->summarize(read, mate);
    }

    list<GraphAlignment> readAlignments = align(read);
    list<GraphAlignment> mateAlignments = align(mate);

    int numMatchingBases = static_cast<int>(read.sequence().length() / 7.5);
    numMatchingBases = std::max(numMatchingBases, 10);
    LinearAlignmentParameters parameters;
    const int kMinNonRepeatAlignmentScore = numMatchingBases * parameters.matchScore;

    if (!checkIfComesFromGraphLocus(readAlignments, mateAlignments, kMinNonRepeatAlignmentScore))
    {
        if (offtargetProcessor_ != nullptr)
        {
            offtargetProcessor_->summarize(read, mate);
        }

        return;
    }

    if (!readAlignments.empty() && !mateAlignments.empty())
    {
        for (auto& feature : features_)
        {
            feature->summarize(read, readAlignments, mate, mateAlignments);
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

void GraphModel::addGraphFeature(GraphFeature* feature) { features_.push_back(feature); }

std::vector<Feature*> GraphModel::modelFeatures()
{
    std::vector<Feature*> modelFeatures;
    for (const auto& feature : features_)
    {
        modelFeatures.push_back(feature);
    }

    if (offtargetProcessor_ != nullptr)
    {
        modelFeatures.push_back(offtargetProcessor_);
    }

    return modelFeatures;
}

void GraphModel::addOfftargetReadProcessor(LinearFeature* offtargetProcessor)
{
    if (offtargetProcessor_ != nullptr)
    {
        throw std::runtime_error("Multiple rare repeats at the same locus are not allowed");
    }

    offtargetProcessor_ = offtargetProcessor;
}

}