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

#include "output/VcfWriter.hh"

#include <deque>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/algorithm/string/join.hpp>

#include "output/VcfHeader.hh"
#include "output/VcfWriterHelpers.hh"
#include "stats/ReadSupportCalculator.hh"

using std::deque;
using std::map;
using std::ostream;
using std::pair;
using std::set;
using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

void writeBodyHeader(const string& sampleName, ostream& out)
{
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sampleName << "\n";
}

std::ostream& operator<<(std::ostream& out, VcfWriter& vcfWriter)
{
    outputVcfHeader(vcfWriter.regionCatalog_, vcfWriter.sampleFindings_, out);
    writeBodyHeader(vcfWriter.sampleId_, out);
    vcfWriter.writeBody(out);
    return out;
}

VcfWriter::VcfWriter(
    std::string sampleId, Reference& reference, const LocusCatalog& regionCatalog, const SampleFindings& sampleFindings)
    : sampleId_(std::move(sampleId))
    , reference_(reference)
    , regionCatalog_(regionCatalog)
    , sampleFindings_(sampleFindings)
{
}

void VcfWriter::writeBody(ostream& out)
{
    const std::vector<VcfWriter::LocusIdAndVariantId>& ids = VcfWriter::getSortedIdPairs();

    for (const auto& pair : ids)
    {
        const string& locusId = pair.first;
        auto locusSpecPtr = regionCatalog_.at(locusId);

        const LocusFindings& locusFindings = sampleFindings_.at(locusId);
        const string& variantId = pair.second;
        const auto& variantFindings = locusFindings.findingsForEachVariant.at(variantId);

        assert(locusFindings.optionalStats);
        const double locusDepth = locusFindings.optionalStats->depth();

        GraphVariantVcfWriter variantWriter(reference_, locusSpecPtr, locusDepth, out);
        variantFindings->accept(variantWriter);
    }
}

std::vector<VcfWriter::LocusIdAndVariantId> VcfWriter::getSortedIdPairs()
{
    using VariantTuple = std::tuple<int32_t, int64_t, int64_t, LocusIdAndVariantId>;
    std::vector<VariantTuple> tuples;

    for (const auto& locusIdAndFindings : sampleFindings_)
    {
        const string& locusId = locusIdAndFindings.first;
        auto locusSpec = regionCatalog_.at(locusId);

        const LocusFindings& locusFindings = locusIdAndFindings.second;

        for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
        {
            const string& variantId = variantIdAndFindings.first;
            const GenomicRegion& variantLocation = locusSpec->getVariantLocationById(variantId);
            tuples.emplace_back(
                variantLocation.contigIndex(), variantLocation.start(), variantLocation.end(),
                LocusIdAndVariantId(locusId, variantId));
        }
    }

    std::sort(tuples.begin(), tuples.end());
    std::vector<VcfWriter::LocusIdAndVariantId> idPairs;
    for (const auto& t : tuples)
    {
        idPairs.emplace_back(std::get<3>(t));
    }

    return idPairs;
}

static string createRepeatAlleleSymbol(int repeatSize) { return "<STR" + std::to_string(repeatSize) + ">"; }

static string computeAltSymbol(const RepeatGenotype& genotype, int referenceSizeInUnits)
{
    vector<string> alleleEncodings;

    if (genotype.shortAlleleSizeInUnits() != referenceSizeInUnits)
    {
        alleleEncodings.push_back(createRepeatAlleleSymbol(genotype.shortAlleleSizeInUnits()));
    }

    if (genotype.longAlleleSizeInUnits() != referenceSizeInUnits
        && genotype.shortAlleleSizeInUnits() != genotype.longAlleleSizeInUnits())
    {
        alleleEncodings.push_back(createRepeatAlleleSymbol(genotype.longAlleleSizeInUnits()));
    }

    if (alleleEncodings.empty())
    {
        return ".";
    }

    return boost::algorithm::join(alleleEncodings, ",");
}

static string computeInfoFields(const GraphVariantSpec& variantSpec, const string& repeatUnit)
{
    const auto& referenceLocus = variantSpec.location();
    const int referenceSizeInBp = referenceLocus.length();
    const int referenceSizeInUnits = referenceSizeInBp / repeatUnit.length();

    vector<string> fields;
    fields.push_back("END=" + std::to_string(referenceLocus.end()));
    fields.push_back("REF=" + std::to_string(referenceSizeInUnits));
    fields.push_back("RL=" + std::to_string(referenceSizeInBp));
    fields.push_back("RU=" + repeatUnit);
    fields.push_back("VARID=" + variantSpec.id());
    fields.push_back("REPID=" + variantSpec.id());

    return boost::algorithm::join(fields, ";");
}

static ReadType determineSupportType(const CountTable& spanningCounts, const CountTable& flankingCounts, int repeatSize)
{
    if (spanningCounts.countOf(repeatSize) != 0)
    {
        return ReadType::kSpanning;
    }
    else if (flankingCounts.countOf(repeatSize) != 0)
    {
        return ReadType::kFlanking;
    }

    return ReadType::kRepeat;
}

static string
computeAlleleFields(const GraphVariantSpec& variantSpec, const string& motif, const StrFindings& strFindings)
{
    const auto referenceLocus = variantSpec.location();
    const int referenceSizeInBp = referenceLocus.length();
    const int referenceSizeInUnits = referenceSizeInBp / motif.length();
    const RepeatGenotype& genotype = *strFindings.optionalGenotype();

    ReadSupportCalculator readSupportCalculator(
        strFindings.countsOfSpanningReads(), strFindings.countsOfFlankingReads(), strFindings.countsOfInrepeatReads());

    VcfAlleleFields alleleFields(referenceSizeInUnits);

    const int shortAlleleSize = genotype.shortAlleleSizeInUnits();
    ReadType shortAlleleSupportType = determineSupportType(
        strFindings.countsOfSpanningReads(), strFindings.countsOfFlankingReads(), shortAlleleSize);

    alleleFields.addAlleleInfo(
        shortAlleleSize, shortAlleleSupportType, genotype.shortAlleleSizeInUnitsCi(),
        readSupportCalculator.getCountOfConsistentSpanningReads(shortAlleleSize),
        readSupportCalculator.getCountOfConsistentFlankingReads(shortAlleleSize),
        readSupportCalculator.getCountOfConsistentRepeatReads(shortAlleleSize));

    if (genotype.numAlleles() == 2)
    {
        const int longAlleleSize = genotype.longAlleleSizeInUnits();
        ReadType longAlleleSupportType = determineSupportType(
            strFindings.countsOfSpanningReads(), strFindings.countsOfFlankingReads(), longAlleleSize);

        alleleFields.addAlleleInfo(
            genotype.longAlleleSizeInUnits(), longAlleleSupportType, genotype.longAlleleSizeInUnitsCi(),
            readSupportCalculator.getCountOfConsistentSpanningReads(longAlleleSize),
            readSupportCalculator.getCountOfConsistentFlankingReads(longAlleleSize),
            readSupportCalculator.getCountOfConsistentRepeatReads(longAlleleSize));
    }

    return alleleFields.encode();
}

void GraphVariantVcfWriter::visit(const StrFindings& strFindings)
{
    if (!strFindings.optionalGenotype())
    {
        return;
    }

    auto graphLocusSpecPtr = std::static_pointer_cast<GraphLocusSpec>(locusSpecPtr_);
    const GraphVariantSpec& variantSpec = graphLocusSpecPtr->getVariantById(strFindings.variantId());

    const auto& genotype = *(strFindings.optionalGenotype());

    const auto& referenceLocus = variantSpec.location();
    const auto repeatNodeId = variantSpec.nodes().front();
    const string& repeatUnit = graphLocusSpecPtr->graph().nodeSeq(repeatNodeId);

    const int referenceSizeInUnits = referenceLocus.length() / repeatUnit.length();

    const string altSymbol = computeAltSymbol(genotype, referenceSizeInUnits);
    const string infoFields = computeInfoFields(variantSpec, repeatUnit);
    const string alleleFields = computeAlleleFields(variantSpec, repeatUnit, strFindings);
    const string sampleFields = alleleFields + ":" + std::to_string(locusDepth_);

    const int posPreceedingRepeat1Based = referenceLocus.start();
    const auto& contigName = reference_.contigInfo().getContigName(referenceLocus.contigIndex());
    const string leftFlankingBase
        = reference_.getSequence(contigName, referenceLocus.start() - 1, referenceLocus.start());

    vector<string> vcfRecordElements
        = { contigName, to_string(posPreceedingRepeat1Based),  ".",         leftFlankingBase, altSymbol, ".", "PASS",
            infoFields, "GT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC", sampleFields };

    out_ << boost::algorithm::join(vcfRecordElements, "\t") << std::endl;
}

void GraphVariantVcfWriter::visit(const CnvVariantFindings& cnvFindings)
{
    auto cnvLocusSpecPtr = std::static_pointer_cast<CnvLocusSpec>(locusSpecPtr_);
    const auto& variantSpec = cnvLocusSpecPtr->outputVariant();
    const auto& referenceLocus = *(variantSpec.location);
    const auto& contigName = reference_.contigInfo().getContigName(referenceLocus.contigIndex());
    boost::optional<int> copyNumberChange = cnvFindings.copyNumberChange();
    boost::optional<int> absoluteCopyNumber = cnvFindings.absoluteCopyNumber();
    vector<string> vcfRecordElements;
    
    vector<string> fields;
    fields.push_back("END=" + std::to_string(referenceLocus.end()));
    fields.push_back("CNVLEN=" + std::to_string(referenceLocus.end() - referenceLocus.start()));
    fields.push_back("VARID=" + variantSpec.id);

    std::string callFilter = "LowQ";
    std::string copyNumberChangeCall = ".";
    std::string absoluteCopyNumberCall = ".";
    if (absoluteCopyNumber)
    {
        absoluteCopyNumberCall = to_string(*absoluteCopyNumber);
    }
    if (copyNumberChange)
    {
        callFilter = "PASS";
        copyNumberChangeCall = to_string(*copyNumberChange);
    }
    else if (absoluteCopyNumber)
    {
        callFilter = "BaselineFail";
    }
    vcfRecordElements = { contigName, to_string(referenceLocus.start()),
                            ".",        "N",
                            ".",        ".",
                            callFilter,     boost::algorithm::join(fields, ";"),
                            "CN:CNC",      absoluteCopyNumberCall + ":" + copyNumberChangeCall };
    
    out_ << boost::algorithm::join(vcfRecordElements, "\t") << std::endl;
}

void GraphVariantVcfWriter::visit(const SmallVariantFindings& findings)
{
    auto graphLocusSpecPtr = std::static_pointer_cast<GraphLocusSpec>(locusSpecPtr_);
    const auto& variantSpec = graphLocusSpecPtr->getVariantById(findings.variantId());

    const auto& referenceLocus = variantSpec.location();
    const auto& contigName = reference_.contigInfo().getContigName(referenceLocus.contigIndex());
    string refSequence;
    string altSequence;
    int64_t startPosition = -1;

    if ((variantSpec.classification().subtype == GraphVariantClassification::Subtype::kSwap)
        || (variantSpec.classification().subtype == GraphVariantClassification::Subtype::kSMN))
    {
        assert(variantSpec.optionalRefNode());
        const auto refNode = *variantSpec.optionalRefNode();
        const int refNodeIndex = refNode == variantSpec.nodes().front() ? 0 : 1;
        const int altNodeIndex = refNode == variantSpec.nodes().front() ? 1 : 0;

        const auto refNodeId = variantSpec.nodes()[refNodeIndex];
        const auto altNodeId = variantSpec.nodes()[altNodeIndex];

        refSequence = graphLocusSpecPtr->graph().nodeSeq(refNodeId);
        altSequence = graphLocusSpecPtr->graph().nodeSeq(altNodeId);
        // Conversion from 0-based to 1-based coordinates
        startPosition = referenceLocus.start() + 1;
    }
    else if (variantSpec.classification().subtype == GraphVariantClassification::Subtype::kDeletion)
    {
        const string refFlankingBase
            = reference_.getSequence(contigName, referenceLocus.start() - 1, referenceLocus.start());

        const int refNodeId = variantSpec.nodes().front();
        refSequence = refFlankingBase + graphLocusSpecPtr->graph().nodeSeq(refNodeId);
        altSequence = refFlankingBase;
        // Conversion from 0-based to 1-based coordinates
        startPosition = referenceLocus.start();
    }
    else if (variantSpec.classification().subtype == GraphVariantClassification::Subtype::kInsertion)
    {
        const string refFlankingBase
            = reference_.getSequence(contigName, referenceLocus.start() - 1, referenceLocus.start());

        const int altNodeId = variantSpec.nodes().front();
        refSequence = refFlankingBase;
        altSequence = refFlankingBase + graphLocusSpecPtr->graph().nodeSeq(altNodeId);
        // Conversion from 0-based to 1-based coordinates
        startPosition = referenceLocus.start();
    }
    else
    {
        std::ostringstream encoding;
        encoding << variantSpec.classification().type << "/" << variantSpec.classification().subtype;
        throw std::logic_error("Unable to generate VCF record for " + encoding.str());
    }

    const string infoFields = "VARID=" + variantSpec.id();

    vector<string> sampleFields;
    vector<string> sampleValues;

    auto genotype = findings.optionalGenotype();
    sampleFields.emplace_back("GT");
    sampleValues.push_back(genotype ? streamToString(*genotype) : ".");

    sampleFields.emplace_back("AD");
    std::ostringstream adEncoding;
    adEncoding << findings.numRefReads() << "," << findings.numAltReads();
    sampleValues.push_back(adEncoding.str());

    if (variantSpec.classification().subtype == GraphVariantClassification::Subtype::kSMN)
    {
        string dst;
        switch (findings.refAllelePresenceStatus().status)
        {
        case AlleleStatus::kAbsent:
            dst = "+";
            break;
        case AlleleStatus::kPresent:
            dst = "-";
            break;
        case AlleleStatus::kUncertain:
            dst = "?";
            break;
        }
        sampleFields.emplace_back("DST");
        sampleValues.push_back(dst);
        sampleFields.emplace_back("RPL");
        sampleValues.push_back(streamToString(findings.refAllelePresenceStatus().logLikelihoodRatio));
    }

    sampleFields.emplace_back("LC");
    sampleValues.push_back(std::to_string(locusDepth_));

    const string sampleField = boost::algorithm::join(sampleFields, ":");
    const string sampleValue = boost::algorithm::join(sampleValues, ":");

    vector<string> line{
        contigName, streamToString(startPosition), ".", refSequence, altSequence, ".", "PASS", infoFields, sampleField,
        sampleValue
    };
    out_ << boost::algorithm::join(line, "\t") << std::endl;
}
}
