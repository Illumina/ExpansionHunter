//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "graphcore/Graph.hh"
#include "thirdparty/json/json.hpp"

#include "common/Common.hh"
#include "common/GenomicRegion.hh"
#include "common/Parameters.hh"
#include "common/Reference.hh"
#include "locus_spec/LocusSpec.hh"

namespace ehunter
{

struct AnalysisRegions
{
    // Regions in the reference where we expect relevant reads to align
    std::vector<GenomicRegion> regionsWithReads;
    // Regions where additional relevant reads might be found that require filtering or special considerations
    std::vector<GenomicRegion> offtargetRegionsWithReads;
    std::vector<GenomicRegion> statsRegions;
};

struct GraphVariantClassification
{
    enum class Type
    {
        kRepeat,
        kSmallVariant
    };

    enum class Subtype
    {
        kCommonRepeat,
        kRareRepeat,
        kInsertion,
        kDeletion,
        kSwap,
        kSMN
    };

    GraphVariantClassification(Type type, Subtype subtype)
        : type(type)
        , subtype(subtype)
    {
    }

    bool operator==(const GraphVariantClassification& other) const
    {
        return (other.type == type && other.subtype == subtype);
    }

    Type type;
    Subtype subtype;
};

using NodeLocations = std::unordered_map<graphtools::NodeId, GenomicRegion>;

struct ReferenceGraph
{
    ReferenceGraph(graphtools::Graph graph, NodeLocations nodeLocations)
        : graph(std::move(graph))
        , nodeLocations(std::move(nodeLocations))
    {
    }

    graphtools::Graph graph;
    NodeLocations nodeLocations;
};

class GraphVariantSpec
{
public:
    GraphVariantSpec(
        std::string id, GraphVariantClassification classification, GenomicRegion location,
        std::vector<graphtools::NodeId> nodes, boost::optional<graphtools::NodeId> optionalRefNode)
        : id_(std::move(id))
        , classification_(classification)
        , location_(std::move(location))
        , nodes_(std::move(nodes))
        , optionalRefNode_(optionalRefNode)
    {
        assertConsistency();
    }

    const std::string& id() const { return id_; }
    GraphVariantClassification classification() const { return classification_; }
    const GenomicRegion& location() const { return location_; }
    const std::vector<graphtools::NodeId>& nodes() const { return nodes_; }
    const boost::optional<graphtools::NodeId>& optionalRefNode() const { return optionalRefNode_; }

    bool operator==(const GraphVariantSpec& other) const
    {
        return id_ == other.id_ && classification_ == other.classification() && nodes_ == other.nodes_;
    }

    void assertConsistency() const;

private:
    std::string id_;
    GraphVariantClassification classification_;
    GenomicRegion location_;
    std::vector<graphtools::NodeId> nodes_;
    boost::optional<graphtools::NodeId> optionalRefNode_;
};

class GraphLocusSpec : public LocusSpec
{
public:
    GraphLocusSpec(
        std::string locusId, CopyNumberBySex contigCopyNumber, AnalysisRegions analysisRegions,
        ReferenceGraph referenceGraph, GenotyperParameters genotyperParams)
        : LocusSpec(std::move(locusId), contigCopyNumber)
        , referenceGraph_(std::move(referenceGraph))
        , analysisRegions_(std::move(analysisRegions))
        , genotyperParams_(genotyperParams)
    {
    }
    ~GraphLocusSpec() override = default;

    const AnalysisRegions& analysisRegions() const { return analysisRegions_; }
    const graphtools::Graph& graph() const { return referenceGraph_.graph; }
    const NodeLocations& nodeLocations() const { return referenceGraph_.nodeLocations; }
    const GenotyperParameters& genotyperParams() const { return genotyperParams_; }

    const std::vector<GraphVariantSpec>& variants() const { return variants_; }
    const GraphVariantSpec& getVariantById(const std::string& id) const;
    void addVariant(
        std::string id, GraphVariantClassification classification, GenomicRegion referenceLocus,
        std::vector<graphtools::NodeId> nodes, boost::optional<graphtools::NodeId> optionalRefNode);

private:
    ReferenceGraph referenceGraph_;
    std::vector<GraphVariantSpec> variants_;
    AnalysisRegions analysisRegions_;
    GenotyperParameters genotyperParams_;
};

std::ostream& operator<<(std::ostream& out, GraphVariantClassification::Type type);
std::ostream& operator<<(std::ostream& out, GraphVariantClassification::Subtype subtype);
std::ostream& operator<<(std::ostream& out, GraphVariantClassification classification);
std::ostream& operator<<(std::ostream& out, const GraphVariantSpec& spec);

}
