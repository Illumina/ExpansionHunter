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

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "spdlog/fmt/ostr.h"
#include "spdlog/spdlog.h"

#include "common/Parameters.hh"
#include "common/WorkflowContext.hh"
#include "input/CatalogLoading.hh"
#include "input/ParameterLoading.hh"
#include "input/SampleStats.hh"
#include "output/BamletWriter.hh"
#include "output/JsonWriter.hh"
#include "output/VcfWriter.hh"
//#include "region_analysis/VariantFindings.hh"
#include "sample_analysis/HtsSeekingSampleAnalysis.hh"
#include "sample_analysis/HtsStreamingSampleAnalysis.hh"
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
    case LogLevel::kTrace:
        spdlog::set_level(spdlog::level::trace);
        break;
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
    spdlog::set_pattern("%Y-%m-%dT%H:%M:%S,[%v]");

    try
    {
        spdlog::info("Starting {}", kProgramVersion);

        auto optionalProgramParameters = tryLoadingProgramParameters(argc, argv);
        if (!optionalProgramParameters)
        {
            return 0;
        }

        ProgramParameters& params = *optionalProgramParameters;
        setLogLevel(params.logLevel());

        // Low-level pass through parameters are stored in the context object.
        initializeWorkflowContext(params.heuristics());
        spdlog::info("Workflow parameter object is initialized with {}", WorkflowContext());

        SampleParameters& sampleParams = params.sample();

        spdlog::info("Analyzing sample {}", sampleParams.id());

        const InputPaths& inputPaths = params.inputPaths();

        spdlog::info("Initializing reference {}", inputPaths.reference());
        FastaReference reference(inputPaths.reference(), extractReferenceContigInfo(inputPaths.htsFile()));

        spdlog::info("Loading variant catalog from disk {}", inputPaths.catalog());
        const RegionCatalog regionCatalog = loadLocusCatalogFromDisk(inputPaths.catalog(), reference);

        const OutputPaths& outputPaths = params.outputPaths();

        BamletWriter bamletWriter(outputPaths.bamlet(), reference.contigInfo(), regionCatalog);

        SampleFindings sampleFindings;
        if (params.analysisMode() == AnalysisMode::kSeeking)
        {
            spdlog::info("Running sample analysis in seeking mode");
            sampleFindings = htsSeekingSampleAnalysis(inputPaths, sampleParams.sex(), regionCatalog, bamletWriter);
        }
        else
        {
            spdlog::info("Running sample analysis in streaming mode");
            sampleFindings = htsStreamingSampleAnalysis(inputPaths, sampleParams.sex(), regionCatalog, bamletWriter);
        }

        spdlog::info("Writing output to disk");
        VcfWriter vcfWriter(sampleParams.id(), reference, regionCatalog, sampleFindings);
        writeToFile(outputPaths.vcf(), vcfWriter);

        JsonWriter jsonWriter(sampleParams, reference.contigInfo(), regionCatalog, sampleFindings);
        writeToFile(outputPaths.json(), jsonWriter);
    }
    catch (const std::exception& e)
    {
        spdlog::error(e.what());
        return 1;
    }

    return 0;
}
