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

#pragma once

#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/optional.hpp>

#include "alignment/SoftclippingAligner.hh"
#include "core/Common.hh"
#include "core/GenomicRegion.hh"
#include "core/ReferenceContigInfo.hh"

namespace ehunter
{

enum class AnalysisMode
{
    kSeeking,
    kStreaming
};

enum class LogLevel
{
    kTrace,
    kDebug,
    kInfo,
    kWarn,
    kError
};

class InputPaths
{
public:
    InputPaths(std::string htsFile, std::string reference, std::string catalog)
        : htsFile_(std::move(htsFile))
        , reference_(std::move(reference))
        , catalog_(std::move(catalog))
    {
    }

    const std::string& htsFile() const { return htsFile_; }
    const std::string& reference() const { return reference_; }
    const std::string& catalog() const { return catalog_; }

private:
    std::string htsFile_;
    std::string reference_;
    std::string catalog_;
};

class OutputPaths
{
public:
    OutputPaths(std::string vcf, std::string json, std::string bamlet)
        : vcf_(vcf)
        , json_(json)
        , bamlet_(bamlet)
    {
    }

    const std::string& vcf() const { return vcf_; }
    const std::string& json() const { return json_; }
    const std::string& bamlet() const { return bamlet_; }

private:
    std::string vcf_;
    std::string json_;
    std::string bamlet_;
};

class SampleParameters
{
public:
    SampleParameters(std::string id, Sex sex)
        : id_(std::move(id))
        , sex_(sex)
    {
    }

    const std::string& id() const { return id_; }
    const Sex& sex() const { return sex_; }

private:
    std::string id_;
    Sex sex_;
};

class HeuristicParameters
{
public:
    HeuristicParameters(
        int regionExtensionLength, int minLocusCoverage, int qualityCutoffForGoodBaseCall, bool skipUnaligned,
        const graphtools::AlignerType alignerType, int kmerLenForAlignment = 14, int paddingLength = 10,
        int seedAffixTrimLength = 14, int orientationPredictorKmerLen = 10, int orientationPredictorMinKmerCount = 3)
        : regionExtensionLength_(regionExtensionLength)
        , minLocusCoverage_(minLocusCoverage)
        , qualityCutoffForGoodBaseCall_(qualityCutoffForGoodBaseCall)
        , skipUnaligned_(skipUnaligned)
        , alignerType_(std::move(alignerType))
        , kmerLenForAlignment_(kmerLenForAlignment)
        , paddingLength_(paddingLength)
        , seedAffixTrimLength_(seedAffixTrimLength)
        , orientationPredictorKmerLen_(orientationPredictorKmerLen)
        , orientationPredictorMinKmerCount_(orientationPredictorMinKmerCount)
    {
    }

    int regionExtensionLength() const { return regionExtensionLength_; }
    int minLocusCoverage() const { return minLocusCoverage_; }
    int qualityCutoffForGoodBaseCall() const { return qualityCutoffForGoodBaseCall_; }
    bool skipUnaligned() const { return skipUnaligned_; }
    graphtools::AlignerType alignerType() const { return alignerType_; }
    int kmerLenForAlignment() const { return kmerLenForAlignment_; }
    int paddingLength() const { return paddingLength_; }
    int seedAffixTrimLength() const { return seedAffixTrimLength_; }
    int orientationPredictorKmerLen() const { return orientationPredictorKmerLen_; }
    int orientationPredictorMinKmerCount() const { return orientationPredictorMinKmerCount_; }

private:
    int regionExtensionLength_;
    int minLocusCoverage_;
    int qualityCutoffForGoodBaseCall_;
    bool skipUnaligned_;
    graphtools::AlignerType alignerType_;
    int kmerLenForAlignment_;
    int paddingLength_;
    int seedAffixTrimLength_;
    int orientationPredictorKmerLen_;
    int orientationPredictorMinKmerCount_;
};

// Per-locus parameters (settable from variant catalog) controlling genotyping
struct GenotyperParameters
{
    GenotyperParameters(int minLocusCoverage)
        : minLocusCoverage(minLocusCoverage)
    {
    }
    // Base error rate assumed in SNV key-allele genotyping model
    double errorRate = 0.02;
    // Threshold to call SNV key-allele confidently present / absent
    double likelihoodRatioThreshold = 10000;
    // Minimal estimated locus coverage (depth) to attempt genotyping
    double minLocusCoverage;
    // Minimal number of reads spanning a variant breakpoint
    int minBreakpointSpanningReads = 5;
};

class ProgramParameters
{
public:
    ProgramParameters(
        InputPaths inputPaths, OutputPaths outputPaths, SampleParameters sample, HeuristicParameters heuristics,
        AnalysisMode analysisMode, LogLevel logLevel, const int initThreadCount, const bool initDisableBamletOutput)
        : threadCount(initThreadCount)
        , disableBamletOutput(initDisableBamletOutput)
        , inputPaths_(std::move(inputPaths))
        , outputPaths_(std::move(outputPaths))
        , sample_(std::move(sample))
        , heuristics_(std::move(heuristics))
        , analysisMode_(analysisMode)
        , logLevel_(logLevel)
    {
    }

    const InputPaths& inputPaths() const { return inputPaths_; }
    const OutputPaths& outputPaths() const { return outputPaths_; }
    const SampleParameters& sample() const { return sample_; }
    const HeuristicParameters& heuristics() const { return heuristics_; }
    AnalysisMode analysisMode() const { return analysisMode_; }
    LogLevel logLevel() const { return logLevel_; }

    int threadCount;
    bool disableBamletOutput;

private:
    InputPaths inputPaths_;
    OutputPaths outputPaths_;
    SampleParameters sample_;
    HeuristicParameters heuristics_;
    AnalysisMode analysisMode_;
    LogLevel logLevel_;
};

}
