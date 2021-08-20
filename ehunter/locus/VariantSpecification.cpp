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

#include "locus/VariantSpecification.hh"

#include <sstream>

using std::string;
using std::to_string;

namespace ehunter
{

void VariantSpecification::assertConsistency() const
{
    const bool variantIsRepeat = classification_.type == VariantType::kRepeat;
    const bool variantIsDeletionOrSwap = classification_.type == VariantType::kSmallVariant
        && (classification_.subtype == VariantSubtype::kDeletion || classification_.subtype == VariantSubtype::kSwap
            || classification_.subtype == VariantSubtype::kSMN);
    const bool variantIsInsertion
        = classification_.type == VariantType::kSmallVariant && classification_.subtype == VariantSubtype::kInsertion;

    bool variantIsValid = false;

    if (variantIsRepeat)
    {
        variantIsValid = classification_.subtype == VariantSubtype::kCommonRepeat
            || classification_.subtype == VariantSubtype::kRareRepeat;
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

std::ostream& operator<<(std::ostream& out, VariantType type)
{
    switch (type)
    {
    case VariantType::kSmallVariant:
        out << "SmallVariant";
        break;
    case VariantType::kRepeat:
        out << "Repeat";
        break;
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, VariantSubtype subtype)
{
    switch (subtype)
    {
    case VariantSubtype::kRareRepeat:
        out << "RareRepeat";
        break;
    case VariantSubtype::kCommonRepeat:
        out << "Repeat";
        break;
    case VariantSubtype::kDeletion:
        out << "Deletion";
        break;
    case VariantSubtype::kInsertion:
        out << "Insertion";
        break;
    case VariantSubtype::kSwap:
        out << "Swap";
        break;
    case VariantSubtype::kSMN:
        out << "SMN";
        break;
    }

    return out;
}

std::ostream& operator<<(std::ostream& out, VariantClassification classification)
{
    out << classification.type << "/" << classification.subtype;
    return out;
}

std::ostream& operator<<(std::ostream& out, const VariantSpecification& variantSpec)
{
    const string refNodeEncoding = variantSpec.optionalRefNode() ? to_string(*variantSpec.optionalRefNode()) : "None";
    out << "ID=" << variantSpec.id() << ";classification=" << variantSpec.classification()
        << ";ReferenceLocus=" << variantSpec.referenceLocus() << ";optionalRefNode=" << refNodeEncoding;

    return out;
}

}
