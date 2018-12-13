//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#pragma once

// cppcheck-suppress missingInclude
#include "nlohmann/json.hpp"

#include "graphcore/Graph.hh"
#include "graphcore/GraphReferenceMapping.hh"

using namespace graphtools;

namespace graphIO
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
