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

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

#include "reads/Read.hh"
#include "workflow/ModelFeature.hh"

namespace ehunter
{

class GraphModel;

class GraphFeature : public ModelFeature
{
public:
    GraphFeature(std::shared_ptr<GraphModel> modelPtr, std::vector<graphtools::NodeId> nodeIds);
    ~GraphFeature() override = default;

    using Alignments = std::list<graphtools::GraphAlignment>;
    virtual void process(const Read& read, const Alignments& readAligns, const Read& mate, const Alignments& mateAligns)
        = 0;

    std::shared_ptr<RegionModel> model() override;

protected:
    std::shared_ptr<GraphModel> modelPtr_;
    std::vector<graphtools::NodeId> nodeIds_;
};

}
