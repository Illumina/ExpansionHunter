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

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "graphcore/Graph.hh"

#include "core/Common.hh"
#include "core/GenomicRegion.hh"
#include "core/Parameters.hh"

namespace ehunter
{

enum class VariantType
{
    kRepeat,
    kSmallVariant
};

enum class VariantSubtype
{
    kCommonRepeat,
    kRareRepeat,
    kInsertion,
    kDeletion,
    kSwap,
    kSMN
};

struct VariantClassification
{
    VariantClassification(VariantType type, VariantSubtype subtype)
        : type(type)
        , subtype(subtype)
    {
    }

    bool operator==(const VariantClassification& other) const
    {
        return (other.type == type && other.subtype == subtype);
    }

    VariantType type;
    VariantSubtype subtype;
};

class VariantSpecification
{
public:
    VariantSpecification(
        std::string id, VariantClassification classification, GenomicRegion referenceLocus,
        std::vector<graphtools::NodeId> nodes, boost::optional<graphtools::NodeId> optionalRefNode)
        : id_(std::move(id))
        , classification_(classification)
        , referenceLocus_(std::move(referenceLocus))
        , nodes_(std::move(nodes))
        , optionalRefNode_(optionalRefNode)
    {
        assertConsistency();
    }

    const std::string& id() const { return id_; }
    VariantClassification classification() const { return classification_; }
    const GenomicRegion& referenceLocus() const { return referenceLocus_; }
    const std::vector<graphtools::NodeId>& nodes() const { return nodes_; }
    const boost::optional<graphtools::NodeId>& optionalRefNode() const { return optionalRefNode_; }

    bool operator==(const VariantSpecification& other) const
    {
        return id_ == other.id_ && classification_ == other.classification() && nodes_ == other.nodes_;
    }

    void assertConsistency() const;

private:
    std::string id_;
    VariantClassification classification_;
    GenomicRegion referenceLocus_;
    std::vector<graphtools::NodeId> nodes_;
    boost::optional<graphtools::NodeId> optionalRefNode_;
};

std::ostream& operator<<(std::ostream& out, VariantType type);
std::ostream& operator<<(std::ostream& out, VariantSubtype subtype);
std::ostream& operator<<(std::ostream& out, VariantClassification classification);
std::ostream& operator<<(std::ostream& out, const VariantSpecification& variantSpec);

}
