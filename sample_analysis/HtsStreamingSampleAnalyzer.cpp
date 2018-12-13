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

#include "sample_analysis/HtsStreamingSampleAnalyzer.hh"

#include <memory>
#include <unordered_map>

#include "region_analysis/RegionAnalyzer.hh"
#include "sample_analysis/HtsFileStreamer.hh"
#include "sample_analysis/HtsHelpers.hh"
#include "sample_analysis/LocationBasedDispatcher.hh"

using std::map;
using std::string;
using std::vector;

using std::unordered_map;

namespace ehunter
{

SampleFindings htslibStreamingSampleAnalyzer(
    const InputPaths& inputPaths, const SampleParameters& sampleParams, const HeuristicParameters& heuristicParams,
    const RegionCatalog& regionCatalog, std::ostream& alignmentStream)
{
    vector<std::unique_ptr<RegionAnalyzer>> locusAnalyzers
        = initializeRegionAnalyzers(regionCatalog, sampleParams, heuristicParams, alignmentStream);
    LocationBasedDispatcher locationBasedDispatcher(locusAnalyzers, heuristicParams.regionExtensionLength());

    htshelpers::HtsFileStreamer readStreamer(inputPaths.htsFile());
    while (readStreamer.trySeekingToNextPrimaryAlignment() && readStreamer.isStreamingAlignedReads())
    {
        locationBasedDispatcher.dispatch(
            readStreamer.currentReadChrom(), readStreamer.currentReadPosition(), readStreamer.currentMateChrom(),
            readStreamer.currentMatePosition(), readStreamer.decodeRead());
    }

    SampleFindings sampleFindings;
    for (auto& locusAnalyzer : locusAnalyzers)
    {
        auto locusFindings = locusAnalyzer->genotype();
        sampleFindings.emplace(std::make_pair(locusAnalyzer->regionId(), std::move(locusFindings)));
    }

    return sampleFindings;
}

}
