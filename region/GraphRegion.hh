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
#include <vector>

#include "region/Region.hh"

#include "graphalign/GappedAligner.hh"
#include "graphcore/Graph.hh"

#include "common/Parameters.hh"
#include "filtering/OrientationPredictor.hh"

namespace ehunter
{

class GraphFeature;

class GraphRegion : public Region
{
public:
    using SPtr = std::shared_ptr<GraphRegion>;
    explicit GraphRegion(std::string locusId, graphtools::Graph graph, const HeuristicParameters& heuristics);
    ~GraphRegion() override = default;

    void analyze(Read read, boost::optional<Read> mate) override;
    const graphtools::Graph& graph() const { return graph_; }

private:
    std::list<graphtools::GraphAlignment> align(Read& read) const;

    graphtools::Graph graph_;
    graphtools::GappedGraphAligner aligner_;
    OrientationPredictor orientationPredictor_;

    std::vector<GraphFeature*> features_;
};

class GraphFeature
{
public:
    using SPtr = std::shared_ptr<GraphFeature>;

    explicit GraphFeature(GraphRegion::SPtr regionModelPtr, std::vector<graphtools::NodeId> nodeIds)
        : regionModelPtr_(std::move(regionModelPtr))
        , nodeIds_(std::move(nodeIds))
    {
    }
    virtual ~GraphFeature() = default;

    using Alignments = std::list<graphtools::GraphAlignment>;
    virtual void process(const Read& read, const Alignments& readAligns, const Read& mate, const Alignments& mateAligns)
        = 0;
    GraphRegion::SPtr regionModelPtr() const { return regionModelPtr_; }

protected:
    GraphRegion::SPtr regionModelPtr_;
    std::vector<graphtools::NodeId> nodeIds_;
};

}
