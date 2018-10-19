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
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "common/common.h"
#include "common/genomic_region.h"

class RegionBlueprintComponent
{
public:
    enum class Type
    {
        kFlank,
        kRepeat,
        kInterruption
    };

    enum class Rarity
    {
        kCommon,
        kRare
    };

    RegionBlueprintComponent(const std::string id, const std::string& sequence, Type type, Rarity rarity)
        : id_(id)
        , sequence_(sequence)
        , type_(type)
        , rarity_(rarity)
    {
    }

    const std::string& id() const { return id_; }
    const std::string& sequence() const { return sequence_; }
    Type type() const { return type_; }
    Rarity rarity() const { return rarity_; }

    void setReferenceRegion(const Region& region) { referenceRegion_ = region; }
    boost::optional<Region> referenceRegion() const { return referenceRegion_; }

    bool operator==(const RegionBlueprintComponent& other) const
    {
        return id_ == other.id_ && sequence_ == other.sequence_ && type_ == other.type_
            && referenceRegion_ == other.referenceRegion_;
    }

private:
    std::string id_;
    std::string sequence_;
    Type type_;
    Rarity rarity_;
    boost::optional<Region> referenceRegion_;
};

class RegionBlueprint
{
public:
    using size_type = size_t;
    using const_iterator = std::vector<RegionBlueprintComponent>::const_iterator;
    const_iterator begin() const { return components_.begin(); }
    const_iterator end() const { return components_.end(); }
    const RegionBlueprintComponent& front() const { return components_.front(); }
    const RegionBlueprintComponent& back() const { return components_.back(); }
    size_type size() const { return components_.size(); }

    RegionBlueprint(
        const std::string& leftFlank, const std::string& regionStructureEncoding, const std::string& rightFlank,
        const std::vector<std::string>& repeatIds, const std::vector<Region>& repeatReferenceRegions,
        const std::vector<RegionBlueprintComponent::Rarity>& repeatRarities);

    bool operator==(const RegionBlueprint& other) const { return components_ == other.components_; }

private:
    std::vector<RegionBlueprintComponent> components_;
};

std::ostream& operator<<(std::ostream& out, RegionBlueprintComponent::Type componentType);
std::ostream& operator<<(std::ostream& out, RegionBlueprintComponent::Rarity componentRarity);
std::ostream& operator<<(std::ostream& out, const RegionBlueprintComponent& component);

std::vector<LabeledSequence<RegionBlueprintComponent::Type>> decodeRegionBlueprintSequence(const std::string& encoding);
