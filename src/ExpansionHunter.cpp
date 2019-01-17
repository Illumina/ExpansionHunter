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

#include <fstream>
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
#include "output/BamletWriter.hh"
#include "output/JsonWriter.hh"
#include "output/VcfWriter.hh"
#include "region_analysis/VariantFindings.hh"
#include "sample_analysis/HtsSeekingSampleAnalyzer.hh"
#include "sample_analysis/HtsStreamingSampleAnalyzer.hh"
#include "src/Version.hh"

namespace spd = spdlog;

using namespace ehunter;

template <typename T> static void writeToFile(std::string fileName, T streamable)
{
    std::ofstream out;

    out.open(fileName.c_str());
    if (!out.is_open())
    {
        throw std::runtime_error("Failed to open " + fileName + " for writing (" + strerror(errno) + ")");
    }

    out << streamable;
}

void setLogLevel(LogLevel logLevel)
{
    switch (logLevel)
    {
    case LogLevel::kDebug:
        spdlog::set_level(spdlog::level::debug);
        break;
    case LogLevel::kInfo:
        spdlog::set_level(spdlog::level::info);
        break;
    case LogLevel::kWarn:
        spdlog::set_level(spdlog::level::warn);
        break;
    case LogLevel::kError:
        spdlog::set_level(spdlog::level::err);
        break;
    }
}

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

        setLogLevel(params.logLevel());

        SampleParameters& sampleParams = params.sample();

        console->info("Analyzing sample {}", sampleParams.id());
        console->info("Read length is set to {}", sampleParams.readLength());

        const InputPaths& inputPaths = params.inputPaths();

        console->info("Initializing reference {}", inputPaths.reference());
        FastaReference reference(inputPaths.reference(), extractReferenceContigInfo(inputPaths.htsFile()));

        console->info("Loading variant catalog from disk {}", inputPaths.catalog());
        const HeuristicParameters& heuristicParams = params.heuristics();
        const RegionCatalog regionCatalog
            = loadLocusCatalogFromDisk(inputPaths.catalog(), sampleParams.sex(), heuristicParams, reference);

        console->info("Running sample analysis");
        const OutputPaths& outputPaths = params.outputPaths();

        BamletWriter bamletWriter(outputPaths.bamlet(), reference.contigInfo(), regionCatalog);

        SampleFindings sampleFindings;
        if (isBamFile(inputPaths.htsFile()))
        {
            sampleFindings
                = htsSeekingSampleAnalysis(inputPaths, sampleParams, heuristicParams, regionCatalog, bamletWriter);
        }
        else
        {
            sampleFindings
                = htslibStreamingSampleAnalyzer(inputPaths, sampleParams, heuristicParams, regionCatalog, bamletWriter);
        }

        console->info("Writing output to disk");
        VcfWriter vcfWriter(sampleParams, reference, regionCatalog, sampleFindings);
        writeToFile(outputPaths.vcf(), vcfWriter);

        JsonWriter jsonWriter(sampleParams, reference.contigInfo(), regionCatalog, sampleFindings);
        writeToFile(outputPaths.json(), jsonWriter);
    }
    catch (const std::exception& e)
    {
        console->error(e.what());
        return 1;
    }

    return 0;
}
