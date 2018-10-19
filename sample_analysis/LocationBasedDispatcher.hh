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

#include <memory>
#include <string>

#include "reads/read.h"
#include "region_analysis/RegionAnalyzer.hh"
#include "sample_analysis/LocationBasedAnalyzerFinder.hh"

class LocationBasedDispatcher
{
public:
    LocationBasedDispatcher(std::vector<std::unique_ptr<RegionAnalyzer>>& locusAnalyzers, int searchRadius);
    void dispatch(
        const std::string& readChrom, int32_t readPosition, const std::string& mateChrom, int32_t matePosition,
        reads::Read read);

private:
    LocationBasedAnalyzerFinder locationBasedAnalyzerFinder_;

    using ReadCatalog = std::unordered_map<std::string, reads::Read>;
    ReadCatalog unpairedReads_;
};
