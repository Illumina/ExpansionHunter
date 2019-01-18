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

#include "common/HtsHelpers.hh"
#include "region_analysis/RegionAnalyzer.hh"
#include "sample_analysis/HtsFileStreamer.hh"
#include "sample_analysis/LocationBasedDispatcher.hh"

using graphtools::AlignmentWriter;
using std::map;
using std::string;
using std::unordered_map;
using std::vector;

namespace ehunter
{

SampleFindings htslibStreamingSampleAnalyzer(
    const InputPaths& inputPaths, const SampleParameters& sampleParams, const HeuristicParameters& heuristicParams,
    const RegionCatalog& regionCatalog, AlignmentWriter& bamletWriter)
{
    vector<std::unique_ptr<RegionAnalyzer>> locusAnalyzers
        = initializeRegionAnalyzers(regionCatalog, heuristicParams, bamletWriter);
    LocationBasedDispatcher locationBasedDispatcher(locusAnalyzers);

    htshelpers::HtsFileStreamer readStreamer(inputPaths.htsFile());
    while (readStreamer.trySeekingToNextPrimaryAlignment() && readStreamer.isStreamingAlignedReads())
    {
        locationBasedDispatcher.dispatch(
            readStreamer.currentReadContigId(), readStreamer.currentReadPosition(), readStreamer.currentMateContigId(),
            readStreamer.currentMatePosition(), readStreamer.decodeRead());
    }

    SampleFindings sampleFindings;
    for (auto& locusAnalyzer : locusAnalyzers)
    {
        auto locusFindings = locusAnalyzer->analyze(sampleParams);
        sampleFindings.emplace(std::make_pair(locusAnalyzer->regionId(), std::move(locusFindings)));
    }

    return sampleFindings;
}

}
