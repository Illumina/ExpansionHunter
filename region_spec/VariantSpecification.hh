//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "graphcore/Graph.hh"

#include "common/Common.hh"
#include "common/GenomicRegion.hh"

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
