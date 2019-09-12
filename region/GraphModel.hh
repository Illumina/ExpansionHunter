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

#include "region/RegionModel.hh"

#include "graphalign/GappedAligner.hh"
#include "graphcore/Graph.hh"

#include "common/GenomicRegion.hh"
#include "common/Parameters.hh"
#include "filtering/OrientationPredictor.hh"
#include "region/GraphFeature.hh"
#include "region/RegionModel.hh"

namespace ehunter
{

class GraphModel : public RegionModel
{
public:
    using SPtr = std::shared_ptr<GraphModel>;
    explicit GraphModel(GenomicRegion referenceRegion, graphtools::Graph graph, const HeuristicParameters& heuristics);
    ~GraphModel() override = default;

    void analyze(Read read, boost::optional<Read> mate) override;
    const graphtools::Graph& graph() const { return graph_; }
    const GenomicRegion& referenceRegion() const { return referenceRegion_; }
    
private:
    std::list<graphtools::GraphAlignment> align(Read& read) const;

    GenomicRegion referenceRegion_;
    graphtools::Graph graph_;
    graphtools::GappedGraphAligner aligner_;
    OrientationPredictor orientationPredictor_;
};

}
