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

#include "locus_spec/GraphLocusSpec.hh"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "spdlog/spdlog.h"
#include "thirdparty/json/json.hpp"

#include "common/Common.hh"
#include "common/Reference.hh"

using boost::optional;
using graphtools::NodeId;
using std::map;
using std::ostream;
using std::string;
using std::to_string;
using std::vector;

using Json = nlohmann::json;

namespace spd = spdlog;

namespace ehunter
{

void GraphLocusSpec::addVariant(
    string id, GraphVariantClassification classification, GenomicRegion referenceLocus, vector<NodeId> nodes,
    optional<NodeId> refNode)
{
    variants_.emplace_back(std::move(id), classification, std::move(referenceLocus), std::move(nodes), refNode);
}

const GraphVariantSpec& GraphLocusSpec::getVariantById(const string& id) const
{
    for (const auto& variant : variants_)
    {
        if (variant.id() == id)
        {
            return variant;
        }
    }

    throw std::logic_error("There is no variant " + id + " in locus " + locusId_);
}

void GraphVariantSpec::assertConsistency() const
{
    const bool variantIsRepeat = classification_.type == GraphVariantClassification::Type::kRepeat;
    const bool variantIsDeletionOrSwap = classification_.type == GraphVariantClassification::Type::kSmallVariant
        && (classification_.subtype == GraphVariantClassification::Subtype::kDeletion
            || classification_.subtype == GraphVariantClassification::Subtype::kSwap
            || classification_.subtype == GraphVariantClassification::Subtype::kSMN);
    const bool variantIsInsertion = classification_.type == GraphVariantClassification::Type::kSmallVariant
        && classification_.subtype == GraphVariantClassification::Subtype::kInsertion;

    bool variantIsValid = false;

    if (variantIsRepeat)
    {
        variantIsValid = classification_.subtype == GraphVariantClassification::Subtype::kCommonRepeat
            || classification_.subtype == GraphVariantClassification::Subtype::kRareRepeat;
    }
    else if (variantIsDeletionOrSwap)
    {
        variantIsValid = optionalRefNode_ != boost::none;
    }
    else if (variantIsInsertion)
    {
        variantIsValid = optionalRefNode_ == boost::none;
    }

    if (!variantIsValid)
    {
        std::ostringstream encoding;
        encoding << *this;
        throw std::logic_error("Definition of variant " + encoding.str() + " is inconsistent");
    }
}

std::ostream& operator<<(std::ostream& out, GraphVariantClassification::Type type)
{
    switch (type)
    {
    case GraphVariantClassification::Type::kSmallVariant:
        out << "SmallVariant";
        break;
    case GraphVariantClassification::Type::kRepeat:
        out << "Repeat";
        break;
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, GraphVariantClassification::Subtype subtype)
{
    switch (subtype)
    {
    case GraphVariantClassification::Subtype::kRareRepeat:
        out << "RareRepeat";
        break;
    case GraphVariantClassification::Subtype::kCommonRepeat:
        out << "Repeat";
        break;
    case GraphVariantClassification::Subtype::kDeletion:
        out << "Deletion";
        break;
    case GraphVariantClassification::Subtype::kInsertion:
        out << "Insertion";
        break;
    case GraphVariantClassification::Subtype::kSwap:
        out << "Swap";
        break;
    case GraphVariantClassification::Subtype::kSMN:
        out << "SMN";
        break;
    }

    return out;
}

std::ostream& operator<<(std::ostream& out, GraphVariantClassification classification)
{
    out << classification.type << "/" << classification.subtype;
    return out;
}

std::ostream& operator<<(std::ostream& out, const GraphVariantSpec& spec)
{
    const string refNodeEncoding = spec.optionalRefNode() ? to_string(*spec.optionalRefNode()) : "None";
    out << "id=" << spec.id() << ";classification=" << spec.classification() << ";ReferenceLocus=" << spec.location()
        << ";optionalRefNode=" << refNodeEncoding;

    return out;
}

}
