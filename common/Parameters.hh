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

#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/optional.hpp>

#include "common/Common.hh"
#include "common/GenomicRegion.hh"

namespace ehunter
{

class Outputs
{
public:
    Outputs(const std::string& vcfPath, const std::string& jsonPath, const std::string& logPath);
    std::ostream& vcf() { return vcf_; }
    std::ostream& json() { return json_; }
    std::ostream& log() { return log_; }

private:
    std::ofstream vcf_;
    std::ofstream json_;
    std::ofstream log_;
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
    OutputPaths(std::string vcf, std::string json, std::string log)
        : vcf_(vcf)
        , json_(json)
        , log_(log)
    {
    }

    const std::string& vcf() const { return vcf_; }
    const std::string& json() const { return json_; }
    const std::string& log() const { return log_; }

private:
    std::string vcf_;
    std::string json_;
    std::string log_;
};

class SampleParameters
{
public:
    SampleParameters(
        std::string id, Sex sex, int readLength,
        boost::optional<double> optionalHaplotypeDepth = boost::optional<double>())
        : id_(std::move(id))
        , sex_(sex)
        , readLength_(readLength)
        , optionalHaplotypeDepth_(optionalHaplotypeDepth)
    {
    }

    const std::string& id() const { return id_; }
    const Sex& sex() const { return sex_; }
    int readLength() const { return readLength_; }

    double haplotypeDepth() const
    {
        if (!optionalHaplotypeDepth_)
        {
            throw std::logic_error("Attempting to access unset depth parameter");
        }
        return *optionalHaplotypeDepth_;
    }

    bool isHaplotypeDepthSet() const { return optionalHaplotypeDepth_.is_initialized(); }
    void setHaplotypeDepth(double haplotypeDepth) { optionalHaplotypeDepth_ = haplotypeDepth; }

private:
    std::string id_;
    Sex sex_;
    int readLength_;
    boost::optional<double> optionalHaplotypeDepth_;
};

class HeuristicParameters
{
public:
    HeuristicParameters(
        bool verboseLogging, int regionExtensionLength, int qualityCutoffForGoodBaseCall, bool skipUnaligned,
        const std::string& alignerType, int kmerLenForAlignment = 14, int paddingLength = 10,
        int seedAffixTrimLength = 5)
        : verboseLogging_(verboseLogging)
        , regionExtensionLength_(regionExtensionLength)
        , qualityCutoffForGoodBaseCall_(qualityCutoffForGoodBaseCall)
        , skipUnaligned_(skipUnaligned)
        , alignerType_(alignerType)
        , kmerLenForAlignment_(kmerLenForAlignment)
        , paddingLength_(paddingLength)
        , seedAffixTrimLength_(seedAffixTrimLength)

    {
    }

    bool verboseLogging() const { return verboseLogging_; }
    int regionExtensionLength() const { return regionExtensionLength_; }
    int qualityCutoffForGoodBaseCall() const { return qualityCutoffForGoodBaseCall_; }
    bool skipUnaligned() const { return skipUnaligned_; }
    const std::string& alignerType() const { return alignerType_; }
    int kmerLenForAlignment() const { return kmerLenForAlignment_; }
    int paddingLength() const { return paddingLength_; }
    int seedAffixTrimLength() const { return seedAffixTrimLength_; }

private:
    bool verboseLogging_;
    int regionExtensionLength_;
    int qualityCutoffForGoodBaseCall_;
    bool skipUnaligned_;
    std::string alignerType_;
    int kmerLenForAlignment_;
    int paddingLength_;
    int seedAffixTrimLength_;
};

class ProgramParameters
{
public:
    ProgramParameters(
        InputPaths inputPaths, OutputPaths outputPaths, SampleParameters sample, HeuristicParameters heuristics)
        : inputPaths_(std::move(inputPaths))
        , outputPaths_(std::move(outputPaths))
        , sample_(std::move(sample))
        , heuristics_(std::move(heuristics))
    {
    }

    const InputPaths& inputPaths() const { return inputPaths_; }
    const OutputPaths& outputPaths() const { return outputPaths_; }
    SampleParameters& sample() { return sample_; }
    const HeuristicParameters& heuristics() const { return heuristics_; }

private:
    InputPaths inputPaths_;
    OutputPaths outputPaths_;
    SampleParameters sample_;
    HeuristicParameters heuristics_;
};

}
