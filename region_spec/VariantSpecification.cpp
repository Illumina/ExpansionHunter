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

#include "region_spec/VariantSpecification.hh"

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
