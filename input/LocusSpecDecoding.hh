//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include <string>
#include <vector>

#include "common/GenomicRegion.hh"
#include "common/Parameters.hh"
#include "common/Reference.hh"
#include "region_spec/LocusSpecification.hh"

namespace ehunter
{

enum class VariantTypeFromUser
{
    kRareRepeat,
    kCommonRepeat,
    kSmallVariant,
    kSMN
};

struct LocusDescriptionFromUser
{
    std::string locusId;
    std::string locusStructure;
    std::vector<GenomicRegion> referenceRegions;
    std::vector<GenomicRegion> targetRegions;
    std::vector<GenomicRegion> offtargetRegions;
    std::vector<VariantTypeFromUser> variantTypesFromUser;
};

void assertValidity(const LocusDescriptionFromUser& userDescription);

LocusSpecification decodeLocusSpecification(
    const LocusDescriptionFromUser& userDescription,
    Sex sampleSex,
    const Reference& reference,
    const HeuristicParameters& heuristicParams);

}
