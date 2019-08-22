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

#include <iostream>
#include <list>
#include <set>
#include <vector>

#include <boost/optional.hpp>

#include "graphalign/GraphAlignment.hh"

#include "classification/AlignmentSummary.hh"
#include "reads/Read.hh"

using graphtools::GraphAlignment;

namespace ehunter
{

class StrAlignmentClassifier
{
public:
    StrAlignmentClassifier(const graphtools::Graph& graph, int repeatNodeId);
    const std::set<graphtools::NodeId>& leftFlankNodeIds() const { return leftFlankNodeIds_; }
    const std::set<graphtools::NodeId>& rightFlankNodeIds() const { return rightFlankNodeIds_; }

    ReadSummaryForStr classifyRead(const std::string& read, const std::list<graphtools::GraphAlignment>& alignments) const;
    boost::optional<StrAlignment> classify(const graphtools::GraphAlignment& alignment) const;

    bool operator==(const StrAlignmentClassifier& other) const;

private:
    bool checkQuality(
        const std::string& read, const graphtools::GraphAlignment& alignment, const StrAlignment& strAlignment) const;

    int repeatNodeId_;
    std::set<graphtools::NodeId> leftFlankNodeIds_;
    std::set<graphtools::NodeId> rightFlankNodeIds_;
};

}
