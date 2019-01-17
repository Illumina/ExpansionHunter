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

#include "input/ParameterLoading.hh"

#include <iostream>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>

#include "input/SampleStats.hh"
#include "src/Version.hh"

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
    optional<int> optionalReadLength;
    optional<double> optionalGenomeCoverage;
    string sampleSexEncoding;

    // Heuristic parameters
    bool verboseLogging;
    string alignerType;
    int regionExtensionLength;
    int qualityCutoffForGoodBaseCall;
    bool skipUnaligned;
};

boost::optional<UserParameters> tryParsingUserParameters(int argc, char** argv)
{
    UserParameters params;

    // clang-format off
    po::options_description usage("Allowed options");
    usage.add_options()
      ("help", "Print help message")
      ("version", "Print version number")
      ("reads", po::value<string>(&params.htsFilePath)->required(), "BAM/CRAM file with aligned reads")
      ("reference", po::value<string>(&params.referencePath)->required(), "FASTA file with reference genome")
      ("variant-catalog", po::value<string>(&params.catalogPath)->required(), "JSON file with variants to genotype")
      ("output-prefix", po::value<string>(&params.outputPrefix)->required(), "Prefix for the output files")
      ("region-extension-length", po::value<int>(&params.regionExtensionLength)->default_value(1000), "How far from on/off-target regions to search for informative reads")
      ("read-length", po::value<int>(), "Read length")
      ("genome-coverage", po::value<double>(), "Read depth on diploid chromosomes")
      ("sex", po::value<string>(&params.sampleSexEncoding)->default_value("female"), "Sex of the sample; must be either male or female")
      ("aligner", po::value<string>(&params.alignerType)->default_value("dag-aligner"), "dag-aligner or path-aligner")
      ("verbose-logging", po::bool_switch(&params.verboseLogging)->default_value(false), "Enable verbose logging");
    // clang-format on

    if (argc == 1)
    {
        std::cerr << usage << std::endl;
        return boost::optional<UserParameters>();
    }

    po::variables_map argumentMap;
    po::store(po::command_line_parser(argc, argv).options(usage).run(), argumentMap);

    if (argumentMap.count("help"))
    {
        std::cerr << usage << std::endl;
        return boost::optional<UserParameters>();
    }

    if (argumentMap.count("version"))
    {
        std::cerr << "Starting " << kProgramVersion << std::endl;
        return boost::optional<UserParameters>();
    }

    po::notify(argumentMap);

    if (argumentMap.count("read-length"))
    {
        params.optionalReadLength = argumentMap["read-length"].as<int>();
    }

    if (argumentMap.count("genome-coverage"))
    {
        params.optionalGenomeCoverage = argumentMap["genome-coverage"].as<double>();
    }

    return params;
}

static void assertWritablePath(const string& pathEncoding)
{
    const fs::path path(pathEncoding);
    const fs::path pathToDirectory = path.parent_path();

    const bool thereIsNoDirectory = pathToDirectory.empty();
    const bool pathLeadsExistingDirectory = fs::is_directory(pathToDirectory);
    const bool filenameIsValid = fs::portable_posix_name(path.filename().string());

    if (!filenameIsValid || (!thereIsNoDirectory && !pathLeadsExistingDirectory))
    {
        throw std::invalid_argument(pathEncoding + " is not a valid output path");
    }
}

static void assertPathToExistingFile(const string& pathEncoding)
{
    const fs::path path(pathEncoding);

    if (fs::is_directory(path) && !fs::exists(path))
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
    // Validate input file paths
    assertPathToExistingFile(userParameters.htsFilePath);
    assertIndexExists(userParameters.htsFilePath);
    assertPathToExistingFile(userParameters.referencePath);
    assertPathToExistingFile(userParameters.catalogPath);

    // Validate output prefix
    assertWritablePath(userParameters.outputPrefix);

    // Validate sample parameters
    if (userParameters.optionalReadLength)
    {
        const int readLength = *userParameters.optionalReadLength;
        const int minAllowedReadLength = 50;
        const int maxAllowedReadLength = 500;
        if (readLength < minAllowedReadLength || readLength > maxAllowedReadLength)
        {
            throw std::invalid_argument(to_string(readLength) + "bp reads are not supported");
        }
    }

    const double kMinDepthAllowed = 10.0;
    if (userParameters.optionalGenomeCoverage && *userParameters.optionalGenomeCoverage < kMinDepthAllowed)
    {
        throw std::invalid_argument("Read depth must be at least " + std::to_string(kMinDepthAllowed));
    }

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
        && userParameters.qualityCutoffForGoodBaseCall > kMaxQualityCutoffForGoodBaseCall)
    {
        const string message = "Base call quality cutoff of " + to_string(userParameters.qualityCutoffForGoodBaseCall)
            + " is not supported; the range of allowed cutoffs is between "
            + to_string(kMinQualityCutoffForGoodBaseCall) + " and " + to_string(kMaxQualityCutoffForGoodBaseCall);
        throw std::invalid_argument(message);
    }
}

SampleParameters decodeSampleParameters(const UserParameters& userParams)
{
    fs::path boostHtsFilePath(userParams.htsFilePath);
    auto sampleId = boostHtsFilePath.stem().string();

    Sex sex = decodeSampleSex(userParams.sampleSexEncoding);

    int readLength
        = userParams.optionalReadLength ? *userParams.optionalReadLength : extractReadLength(userParams.htsFilePath);

    if (userParams.optionalGenomeCoverage)
    {
        const double haplotypeDepth = *userParams.optionalGenomeCoverage / 2;
        return SampleParameters(sampleId, sex, readLength, haplotypeDepth);
    }
    else
    {
        return SampleParameters(sampleId, sex, readLength);
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
    const string logPath = userParams.outputPrefix + ".log";
    OutputPaths outputPaths(vcfPath, jsonPath, logPath);
    SampleParameters sampleParameters = decodeSampleParameters(userParams);
    HeuristicParameters heuristicParameters(
        userParams.verboseLogging, userParams.regionExtensionLength, userParams.qualityCutoffForGoodBaseCall,
        userParams.skipUnaligned, userParams.alignerType);

    return ProgramParameters(inputPaths, outputPaths, sampleParameters, heuristicParameters);
}

}
