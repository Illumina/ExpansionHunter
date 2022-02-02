//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "io/ParameterLoading.hh"

#include <iostream>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>

#include "app/Version.hh"
#include "io/SampleStats.hh"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

using boost::optional;
using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

struct UserParameters
{
    // Input file paths
    string htsFilePath;
    string referencePath;
    string catalogPath;

    // Output prefix
    string outputPrefix;

    // Sample parameters
    string sampleSexEncoding;

    // Heuristic parameters
    string alignerType;
    int regionExtensionLength;
    double minLocusCoverage = 10.0;
    int qualityCutoffForGoodBaseCall = 20;
    bool skipUnaligned;

    string analysisMode;
    string logLevel;
    int threadCount;
    bool disableBamletOutput = false;
};

boost::optional<UserParameters> tryParsingUserParameters(int argc, char** argv)
{
    UserParameters params;

    // clang-format off
    po::options_description basicOptions("Basic options");
    basicOptions.add_options()
        ("help,h", "Print help message")
        ("version,v", "Print version number")
        ("reads", po::value<string>(&params.htsFilePath)->required(), "aligned reads BAM/CRAM file/URL")
        ("reference", po::value<string>(&params.referencePath)->required(), "reference genome FASTA file")
        ("variant-catalog", po::value<string>(&params.catalogPath)->required(), "JSON file with variants to genotype")
        ("output-prefix", po::value<string>(&params.outputPrefix)->required(), "Prefix for the output files")
        ("sex", po::value<string>(&params.sampleSexEncoding)->default_value("female"), "Sex of the sample; must be either male or female")
    ;
    // clang-format on

    // clang-format off
    po::options_description advancedOptions("Advanced options");
    advancedOptions.add_options()
        ("region-extension-length", po::value<int>(&params.regionExtensionLength)->default_value(1000), "How far from on/off-target regions to search for informative reads")
        ("min-locus-coverage", po::value<double>(&params.minLocusCoverage)->default_value(10.0), "Minimum read coverage depth for diploid loci (set to half for loci on haploid chromosomes)")
        ("aligner", po::value<string>(&params.alignerType)->default_value("dag-aligner"), "Graph aligner to use (dag-aligner or path-aligner)")
        ("analysis-mode", po::value<string>(&params.analysisMode)->default_value("seeking"), "Analysis workflow to use (seeking or streaming)")
        ("threads", po::value(&params.threadCount)->default_value(1), "Number of threads to use")
        ("log-level", po::value<string>(&params.logLevel)->default_value("info"), "trace, debug, info, warn, or error")
    ;
    // clang-format on

    // Internal options are added for development and debug activities, these are suppressed from the commandline
    // because they may be removed or changed without deprecation.

    // clang-format off
    po::options_description internalOptions("Internal options (not stable in future releases)");
    internalOptions.add_options()
    ("disable-bamlet-output", "Disable bamlet output")
    ;
    // clang-format on

    po::options_description visibleOptions;
    visibleOptions.add(basicOptions).add(advancedOptions);

    po::options_description cmdlineOptions(visibleOptions);
    cmdlineOptions.add(internalOptions);

    po::variables_map argumentMap;
    po::store(po::command_line_parser(argc, argv).options(cmdlineOptions).run(), argumentMap);

    if ((argc == 1) or argumentMap.count("help"))
    {
        std::cerr << visibleOptions << std::endl;
        return {};
    }

    if (argumentMap.count("version"))
    {
        std::cerr << "Starting " << kProgramVersion << std::endl;
        return {};
    }

    params.disableBamletOutput = argumentMap.count("disable-bamlet-output");

    po::notify(argumentMap);

    return params;
}

static void assertWritablePath(const string& pathEncoding)
{
    const fs::path path(pathEncoding);
    const fs::path pathToDirectory = path.parent_path();

    const bool thereIsNoDirectory = pathToDirectory.empty();
    const bool pathLeadsToExistingDirectory = fs::is_directory(pathToDirectory);
    const bool filenameIsValid = fs::portable_posix_name(path.filename().string());

    if (!filenameIsValid || (!thereIsNoDirectory && !pathLeadsToExistingDirectory))
    {
        throw std::invalid_argument(pathEncoding + " is not a valid output path");
    }
}

static void assertPathToExistingFile(const string& pathEncoding)
{
    const fs::path path(pathEncoding);
    const bool isPathToExistingFile = fs::exists(path) && fs::is_regular_file(path);

    if (!isPathToExistingFile)
    {
        throw std::invalid_argument(pathEncoding + " is not a path to an existing file");
    }
}

static void assertIndexExists(const string& htsFilePath)
{
    const vector<string> kPossibleIndexExtensions = { ".bai", ".csi", ".crai" };

    for (const string& indexExtension : kPossibleIndexExtensions)
    {
        if (fs::exists(htsFilePath + indexExtension))
        {
            return;
        }
    }

    throw std::invalid_argument("Could not find index of " + htsFilePath);
}

void assertValidity(const UserParameters& userParameters)
{
    // Validate analysis Mode:
    if ((userParameters.analysisMode != "seeking") and (userParameters.analysisMode != "streaming"))
    {
        throw std::invalid_argument(userParameters.analysisMode + " is not a valid analysis mode");
    }

    // Validate input file paths
    if (not isURL(userParameters.htsFilePath))
    {
        assertPathToExistingFile(userParameters.htsFilePath);
        if (userParameters.analysisMode != "streaming")
        {
            assertIndexExists(userParameters.htsFilePath);
        }
    }
    assertPathToExistingFile(userParameters.referencePath);
    assertPathToExistingFile(userParameters.catalogPath);

    // Validate output prefix
    assertWritablePath(userParameters.outputPrefix);

    // Validate sample parameters
    if (userParameters.sampleSexEncoding != "female" && userParameters.sampleSexEncoding != "male")
    {
        throw std::invalid_argument(userParameters.sampleSexEncoding + " is not a valid sex encoding");
    }

    // Heuristic parameters
    if (userParameters.alignerType != "dag-aligner" && userParameters.alignerType != "path-aligner")
    {
        throw std::invalid_argument(userParameters.alignerType + " is not a valid aligner type");
    }

    const int kMinExtensionLength = 500;
    const int kMaxExtensionLength = 1500;
    if (userParameters.regionExtensionLength < kMinExtensionLength
        && userParameters.regionExtensionLength > kMaxExtensionLength)
    {
        const string message = "Extension length of size " + to_string(userParameters.regionExtensionLength)
            + " is not supported; the range of allowed extensions is between " + to_string(kMinExtensionLength)
            + " and " + to_string(kMaxExtensionLength);
        throw std::invalid_argument(message);
    }

    const int kMinQualityCutoffForGoodBaseCall = 5;
    const int kMaxQualityCutoffForGoodBaseCall = 40;
    if (userParameters.qualityCutoffForGoodBaseCall < kMinQualityCutoffForGoodBaseCall
        || userParameters.qualityCutoffForGoodBaseCall > kMaxQualityCutoffForGoodBaseCall)
    {
        const string message = "Base call quality cutoff of " + to_string(userParameters.qualityCutoffForGoodBaseCall)
            + " is not supported; the range of allowed cutoffs is between "
            + to_string(kMinQualityCutoffForGoodBaseCall) + " and " + to_string(kMaxQualityCutoffForGoodBaseCall);
        throw std::invalid_argument(message);
    }

    if (userParameters.threadCount < 1)
    {
        const string message = "Thread count cannot be less than 1";
        throw std::invalid_argument(message);
    }
}

SampleParameters decodeSampleParameters(const UserParameters& userParams)
{
    fs::path boostHtsFilePath(userParams.htsFilePath);
    auto sampleId = boostHtsFilePath.stem().string();
    Sex sex = decodeSampleSex(userParams.sampleSexEncoding);
    return SampleParameters(sampleId, sex);
}

AnalysisMode decodeAnalysisMode(const string& encoding)
{
    if (encoding == "streaming")
    {
        return AnalysisMode::kStreaming;
    }
    else if (encoding == "seeking")
    {
        return AnalysisMode::kSeeking;
    }
    else
    {
        throw std::logic_error("Invalid encoding of data input mode '" + encoding + "'");
    }
}

LogLevel decodeLogLevel(const string& encoding)
{
    if (encoding == "trace")
    {
        return LogLevel::kTrace;
    }
    if (encoding == "debug")
    {
        return LogLevel::kDebug;
    }
    else if (encoding == "info")
    {
        return LogLevel::kInfo;
    }
    else if (encoding == "warn")
    {
        return LogLevel::kWarn;
    }
    else if (encoding == "error")
    {
        return LogLevel::kError;
    }
    else
    {
        throw std::logic_error("Invalid encoding of logging level " + encoding);
    }
}

static graphtools::AlignerType decodeAlignerType(const string& alignerType)
{
    if (alignerType == "path-aligner")
    {
        return graphtools::AlignerType::PATH_ALIGNER;
    }
    if (alignerType == "dag-aligner")
    {
        return graphtools::AlignerType::DAG_ALIGNER;
    }
    else
    {
        throw std::logic_error(alignerType + " is not a valid aligner type");
    }
}

boost::optional<ProgramParameters> tryLoadingProgramParameters(int argc, char** argv)
{
    auto optionalUserParameters = tryParsingUserParameters(argc, argv);
    if (!optionalUserParameters)
    {
        return boost::optional<ProgramParameters>();
    }

    const auto& userParams = *optionalUserParameters;
    assertValidity(userParams);

    InputPaths inputPaths(userParams.htsFilePath, userParams.referencePath, userParams.catalogPath);
    const string vcfPath = userParams.outputPrefix + ".vcf";
    const string jsonPath = userParams.outputPrefix + ".json";
    const string bamletPath = userParams.outputPrefix + "_realigned.bam";
    OutputPaths outputPaths(vcfPath, jsonPath, bamletPath);
    SampleParameters sampleParameters = decodeSampleParameters(userParams);
    HeuristicParameters heuristicParameters(
        userParams.regionExtensionLength, userParams.minLocusCoverage, userParams.qualityCutoffForGoodBaseCall,
        userParams.skipUnaligned, decodeAlignerType(userParams.alignerType));

    LogLevel logLevel;
    try
    {
        logLevel = decodeLogLevel(userParams.logLevel);
    }
    catch (std::logic_error&)
    {
        const string message = "Log level must be set to either trace, debug, info, warn, or error";
        throw std::invalid_argument(message);
    }

    AnalysisMode analysisMode;
    try
    {
        analysisMode = decodeAnalysisMode(userParams.analysisMode);
    }
    catch (std::logic_error&)
    {
        const string message = "Analysis mode must be set to either streaming or seeking";
        throw std::invalid_argument(message);
    }

    return ProgramParameters(
        inputPaths, outputPaths, sampleParameters, heuristicParameters, analysisMode, logLevel, userParams.threadCount,
        userParams.disableBamletOutput);
}

}
