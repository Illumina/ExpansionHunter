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

#pragma once

#include <list>
#include <memory>
#include <string>
#include <vector>

#include "graphalign/GappedAligner.hh"
#include "graphcore/Graph.hh"
#include "graphio/AlignmentWriter.hh"

#include "common/GenomicRegion.hh"
#include "common/Parameters.hh"
#include "filtering/OrientationPredictor.hh"
#include "graph_components/ReadClassifier.hh"
#include "workflow/LinearFeature.hh"
#include "workflow/RegionModel.hh"

namespace ehunter
{

class GraphFeature;
class IrrPairDetector;

class GraphModel : public RegionModel
{
public:
    using Alignments = std::list<graphtools::GraphAlignment>;
    using AlignmentWriter = std::shared_ptr<graphtools::AlignmentWriter>;
    using Regions = std::vector<GenomicRegion>;

    enum class Origin
    {
        kTargetRegion,
        kOfftargetRegion,
        kOtherRegion
    };

    struct AlignmentBundle
    {
        AlignmentBundle(Alignments alignments, bool forwardOriented)
            : alignments(std::move(alignments))
            , forwardOriented(forwardOriented)
        {
        }
        Alignments alignments;
        bool forwardOriented;
    };

    GraphModel(
        std::string graphId, const Regions& targetRegions, const Regions& offtargetRegions, graphtools::Graph graph,
        const HeuristicParameters& heuristics, AlignmentWriter alignmentWriter);
    ~GraphModel() override = default;

    void analyze(const MappedRead& read, const MappedRead& mate) override;
    void analyze(MappedRead /*read*/) override {};
    const graphtools::Graph& graph() const { return graph_; }
    void addGraphFeature(GraphFeature* feature);
    void addOfftargetReadProcessor(LinearFeature* offtargetProcessor);
    std::vector<Feature*> modelFeatures() override;

private:
    Origin guessOrigin(const MappedRead& read, const MappedRead& mate);
    Origin guessOrigin(int readLength, const Alignments& readAlignments, const Alignments& mateAlignments);

    AlignmentBundle align(const std::string& sequence) const;
    void analyzeOfftarget(const MappedRead& read, const MappedRead& mate);

    void writeAlignments(
        const MappedRead& read, const AlignmentBundle& readBundle, const MappedRead& mate,
        const AlignmentBundle& mateBundle);

    std::string graphId_;
    Regions targetRegions_;
    std::vector<GraphFeature*> features_;
    LinearFeature* offtargetProcessor_ = nullptr;
    AlignmentWriter alignmentWriter_;

    ReadClassifier readClassifier_;
    graphtools::Graph graph_;
    graphtools::GappedGraphAligner aligner_;
    OrientationPredictor orientationPredictor_;
};

}
