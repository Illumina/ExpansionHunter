//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>
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

#pragma once

// cppcheck-suppress missingInclude
#include "nlohmann/json.hpp"

#include "graphcore/Graph.hh"
#include "graphcore/GraphReferenceMapping.hh"

using namespace graphtools;

namespace graphtools
{

using Json = nlohmann::json;

/**
 * Load JSON from file and parse Graph.
 * Graph can be either directly the top-level object in the json or under 'graph'
 */
Graph loadGraph(std::string const& jsonPath);

/**
 * Create graph from Json representation.
 * @throws if the Json does not represent a valid graph
 */
Graph parseGraph(Json const& jsonGraph);

/**
 * Create Json representation of the graph
 */
Json graphToJson(Graph const& graph);

/**
 * Load Reference mapping from graph description
 * @param jRefmap Json representation of referenceMapping. Must match the graph
 * @param graph Graph to map to reference.
 * @throws if the Json does not represent a valid reference map or does not match the graph
 */
GraphReferenceMapping parseReferenceMapping(Json const& jRefmap, Graph const& graph);
}
