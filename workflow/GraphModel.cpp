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

#include "spdlog/spdlog.h"

#include "alignment/AlignmentFilters.hh"
#include "workflow/GraphFeature.hh"
#include "workflow/IrrPairDetector.hh"

namespace ehunter
{

using graphtools::GappedGraphAligner;
using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::reverseComplement;
using std::list;
using std::string;
using std::vector;

template <typename T> static vector<T> concatenate(vector<T> first, const vector<T>& second)
{
    first.insert(first.end(), second.begin(), second.end());
    return first;
}

using Origin = GraphModel::Origin;
using Alignments = GraphModel::Alignments;
using AlignmentBundle = GraphModel::AlignmentBundle;

GraphModel::GraphModel(
    string graphId, const Regions& targetRegions, const Regions& offtargetRegions, Graph graph,
    const HeuristicParameters& heuristics, AlignmentWriter alignmentWriter)
    : RegionModel(concatenate(targetRegions, offtargetRegions))
    , graphId_(std::move(graphId))
    , targetRegions_(targetRegions)
    , alignmentWriter_(std::move(alignmentWriter))
    , readClassifier_(readExtractionRegions_)
    , graph_(std::move(graph))
    , aligner_(
          &graph_, heuristics.kmerLenForAlignment(), heuristics.paddingLength(), heuristics.seedAffixTrimLength(),
          heuristics.alignerType())
    , orientationPredictor_(&graph_)
{
}

void GraphModel::analyze(const MappedRead& read, const MappedRead& mate)
{
    Origin origin = guessOrigin(read, mate);

    if (origin == Origin::kOfftargetRegion)
    {
        analyzeOfftarget(read, mate);
        return;
    }

    if (origin != Origin::kTargetRegion)
    {
        return;
    }

    auto readBundle = align(read.sequence());
    auto mateBundle = align(mate.sequence());

    const auto& readAlignments = readBundle.alignments;
    const auto& mateAlignments = mateBundle.alignments;

    origin = guessOrigin(read.sequence().length(), readAlignments, mateAlignments);

    if (origin == Origin::kOfftargetRegion)
    {
        analyzeOfftarget(read, mate);
        return;
    }

    if (origin != Origin::kTargetRegion)
    {
        return;
    }

    if (!readAlignments.empty() && !mateAlignments.empty())
    {
        writeAlignments(read, readBundle, mate, mateBundle);

        const auto& readSequence = readBundle.forwardOriented ? read.sequence() : reverseComplement(read.sequence());
        const auto& mateSequence = mateBundle.forwardOriented ? mate.sequence() : reverseComplement(mate.sequence());

        for (auto& feature : features_)
        {
            feature->summarize(readSequence, readAlignments, mateSequence, mateAlignments);
        }
    }
}

GraphModel::Origin GraphModel::guessOrigin(const MappedRead& read, const MappedRead& mate)
{
    switch (readClassifier_.classify(read, mate))
    {
    case RegionProximity::kInside:
        return Origin ::kTargetRegion;
    case RegionProximity::kOverlapsOrNear:
        return Origin::kOtherRegion;
    case RegionProximity::kFar:
        return Origin::kOfftargetRegion;
    }

    // Unreachable, but required to pass linter check
    return Origin::kOtherRegion;
}

Origin GraphModel::guessOrigin(int readLength, const Alignments& readAlignments, const Alignments& mateAlignments)
{
    int numMatchingBases = static_cast<int>(readLength / 7.5);
    numMatchingBases = std::max(numMatchingBases, 10);
    LinearAlignmentParameters parameters;
    const int kMinNonRepeatAlignmentScore = numMatchingBases * parameters.matchScore;

    if (checkIfComesFromGraphLocus(readAlignments, mateAlignments, kMinNonRepeatAlignmentScore))
    {
        return Origin::kTargetRegion;
    }
    else
    {
        return Origin::kOfftargetRegion;
    }
}

void GraphModel::analyzeOfftarget(const MappedRead& read, const MappedRead& mate)
{
    if (offtargetProcessor_ != nullptr)
    {
        offtargetProcessor_->summarize(read, mate);
    }
}

GraphModel::AlignmentBundle GraphModel::align(const string& sequence) const
{
    switch (orientationPredictor_.predict(sequence))
    {
    case OrientationPrediction::kAlignsInOriginalOrientation:
        return AlignmentBundle(aligner_.align(sequence), true);
    case OrientationPrediction::kAlignsInOppositeOrientation:
        return AlignmentBundle(aligner_.align(graphtools::reverseComplement(sequence)), false);
    case OrientationPrediction::kDoesNotAlign:
        return AlignmentBundle({}, true);
    }

    // Unreachable, but required to pass linter check
    return AlignmentBundle({}, true);
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

void GraphModel::writeAlignments(
    const MappedRead& read, const AlignmentBundle& readBundle, const MappedRead& mate,
    const AlignmentBundle& mateBundle)
{
    const auto& readSequence = readBundle.forwardOriented ? read.sequence() : reverseComplement(read.sequence());
    const auto& mateSequence = mateBundle.forwardOriented ? mate.sequence() : reverseComplement(mate.sequence());

    const bool isReadReversed = readBundle.forwardOriented ? read.isReversed() : !read.isReversed();
    const bool isMateReversed = mateBundle.forwardOriented ? mate.isReversed() : !mate.isReversed();

    const auto& readAlignment = readBundle.alignments.front();
    const auto& mateAlignment = mateBundle.alignments.front();

    alignmentWriter_->write(
        graphId_, read.fragmentId(), readSequence, read.isFirstMate(), isReadReversed, isMateReversed, readAlignment);
    alignmentWriter_->write(
        graphId_, mate.fragmentId(), mateSequence, mate.isFirstMate(), isMateReversed, isReadReversed, mateAlignment);
}
}
