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

#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "thirdparty/spdlog/fmt/ostr.h"
#include "thirdparty/spdlog/spdlog.h"

#include "common/Parameters.hh"
#include "input/CatalogLoading.hh"
#include "input/ParameterLoading.hh"
#include "input/SampleStats.hh"
#include "output/JsonWriter.hh"
#include "output/VcfWriter.hh"
#include "region_analysis/VariantFindings.hh"
#include "sample_analysis/HtsSeekingSampleAnalyzer.hh"
#include "sample_analysis/HtsStreamingSampleAnalyzer.hh"
#include "src/Version.hh"

namespace spd = spdlog;

using namespace ehunter;

int main(int argc, char** argv)
{
    auto console = spd::stderr_color_mt("console");
    console->set_pattern("%Y-%m-%dT%H:%M:%S,[%v]");

    try
    {
        console->info("Starting {}", kProgramVersion);

        auto optionalProgramParameters = tryLoadingProgramParameters(argc, argv);
        if (!optionalProgramParameters)
        {
            return 0;
        }
        ProgramParameters& params = *optionalProgramParameters;

        if (params.heuristics().verboseLogging())
        {
            auto verboseConsole = spdlog::stdout_color_mt("verbose");
            verboseConsole->set_pattern("%Y-%m-%dT%H:%M:%S\n%v");
            console->info("Verbose logging enabled");
        }

        SampleParameters& sampleParams = params.sample();

        console->info("Analyzing sample {}", sampleParams.id());
        console->info("Read length is set to {}", sampleParams.readLength());

        const InputPaths& inputPaths = params.inputPaths();

        console->info("Initializing reference {}", inputPaths.reference());
        FastaReference reference(inputPaths.reference());

        console->info("Loading variant catalog from disk {}", inputPaths.catalog());
        const RegionCatalog regionCatalog
            = loadRegionCatalogFromDisk(inputPaths.catalog(), reference, sampleParams.sex());

        console->info("Running sample analysis");
        const HeuristicParameters& heuristicParams = params.heuristics();
        const OutputPaths& outputPaths = params.outputPaths();
        Outputs outputs(outputPaths.vcf(), outputPaths.json(), outputPaths.log());

        SampleFindings sampleFindings;
        if (isBamFile(inputPaths.htsFile()))
        {
            sampleFindings
                = htsSeekingSampleAnalysis(inputPaths, sampleParams, heuristicParams, regionCatalog, outputs.log());
        }
        else
        {
            sampleFindings = htslibStreamingSampleAnalyzer(
                inputPaths, sampleParams, heuristicParams, regionCatalog, outputs.log());
        }

        console->info("Writing output to disk");
        VcfWriter vcfWriter(sampleParams.id(), sampleParams.readLength(), regionCatalog, sampleFindings, reference);
        outputs.vcf() << vcfWriter;

        JsonWriter jsonWriter(sampleParams.id(), sampleParams.readLength(), regionCatalog, sampleFindings);
        outputs.json() << jsonWriter;
    }
    catch (const std::exception& e)
    {
        console->error(e.what());
        return 1;
    }

    return 0;
}
