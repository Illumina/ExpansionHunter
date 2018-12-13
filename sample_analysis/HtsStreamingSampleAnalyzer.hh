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
#include <vector>

#include "common/Parameters.hh"
#include "region_analysis/VariantFindings.hh"
#include "region_spec/LocusSpecification.hh"

namespace ehunter
{

SampleFindings htslibStreamingSampleAnalyzer(
    const InputPaths& inputPaths, const SampleParameters& sampleParams, const HeuristicParameters& heuristicParams,
    const RegionCatalog& regionCatalog, std::ostream& alignmentStream);

}
