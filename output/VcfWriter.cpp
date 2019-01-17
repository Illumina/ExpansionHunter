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

#include "output/VcfWriter.hh"

#include <deque>
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

template <typename T> static string streamToString(const T& streamableObject)
{
    std::stringstream out;
    out << streamableObject;
    return out.str();
}

void writeBodyHeader(const string& sampleName, ostream& out)
{
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sampleName << "\n";
}

std::ostream& operator<<(std::ostream& out, VcfWriter& vcfWriter)
{
    outputVcfHeader(vcfWriter.regionCatalog_, vcfWriter.sampleFindings_, out);
    writeBodyHeader(vcfWriter.sampleParams_.id(), out);
    vcfWriter.writeBody(out);
    return out;
}

VcfWriter::VcfWriter(
    const SampleParameters& sampleParams, Reference& reference, const RegionCatalog& regionCatalog,
    const SampleFindings& sampleFindings)
    : sampleParams_(sampleParams)
    , reference_(reference)
    , regionCatalog_(regionCatalog)
    , sampleFindings_(sampleFindings)
{
}

void VcfWriter::writeBody(ostream& out)
{
    const std::vector<VcfWriter::RegionIdAndVariantId>& ids = VcfWriter::getSortedIdPairs();

    for (const auto& pair : ids)
    {
        const string& regionId = pair.first;
        const LocusSpecification& regionSpec = regionCatalog_.at(regionId);
        const RegionFindings& regionFindings = sampleFindings_.at(regionId);
        
        const string& variantId = pair.second;
        const VariantSpecification& variantSpec = regionSpec.getVariantSpecById(variantId);
        const auto& variantFindings = regionFindings.at(variantId);

        VariantVcfWriter variantWriter(sampleParams_, reference_, regionSpec, variantSpec, out);
        variantFindings->accept(&variantWriter);
    }
}

const std::vector<VcfWriter::RegionIdAndVariantId> VcfWriter::getSortedIdPairs()
{
    using VariantTuple = std::tuple<int32_t, int64_t, int64_t, RegionIdAndVariantId>;
    std::vector<VariantTuple> tuples;

    for (const auto& regionIdAndFindings : sampleFindings_)
    {
        const string& regionId = regionIdAndFindings.first;
        const LocusSpecification& regionSpec = regionCatalog_.at(regionId);
        const RegionFindings& regionFindings = regionIdAndFindings.second;

        for (const auto& variantIdAndFindings : regionFindings)
        {
            const string& variantId = variantIdAndFindings.first;
            const VariantSpecification& variantSpec = regionSpec.getVariantSpecById(variantId);
            tuples.emplace_back(variantSpec.referenceLocus().contigIndex(), variantSpec.referenceLocus().start(), variantSpec.referenceLocus().end(), RegionIdAndVariantId(regionId, variantId));
        }
    }

    std::sort(tuples.begin(), tuples.end());
    std::vector<VcfWriter::RegionIdAndVariantId> idPairs;
    for (const auto& t : tuples)
        idPairs.emplace_back(std::get<3>(t));

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

static string computeInfoFields(const VariantSpecification& variantSpec, const string& repeatUnit)
{
    const auto& referenceLocus = variantSpec.referenceLocus();
    const int referenceSizeInBp = referenceLocus.length();
    const int referenceSizeInUnits = referenceSizeInBp / repeatUnit.length();

    vector<string> fields;
    fields.push_back("END=" + std::to_string(referenceLocus.end()));
    fields.push_back("REF=" + std::to_string(referenceSizeInUnits));
    fields.push_back("RL=" + std::to_string(referenceSizeInUnits));
    fields.push_back("RU=" + repeatUnit);
    fields.push_back("REPID=" + variantSpec.id());

    return boost::algorithm::join(fields, ";");
}

static ReadType
deterimineSupportType(int maxUnitsInRead, CountTable spanningCounts, CountTable flankingCounts, int repeatSize)
{
    if (spanningCounts.countOf(repeatSize) != 0)
    {
        return ReadType::kSpanning;
    }
    else if (flankingCounts.countOf(maxUnitsInRead) != 0)
    {
        return ReadType::kRepeat;
    }

    return ReadType::kFlanking;
}

static string computeSampleFields(
    int readLength, const VariantSpecification& variantSpec, const string& repeatUnit,
    const RepeatFindings& repeatFindings)
{
    const auto referenceLocus = variantSpec.referenceLocus();
    const int referenceSizeInBp = referenceLocus.length();
    const int referenceSizeInUnits = referenceSizeInBp / repeatUnit.length();
    const RepeatGenotype& genotype = *repeatFindings.optionalGenotype();

    const int maxUnitsInRead = std::ceil(readLength / static_cast<double>(repeatUnit.length()));
    ReadSupportCalculator readSupportCalculator(
        maxUnitsInRead, repeatFindings.countsOfSpanningReads(), repeatFindings.countsOfFlankingReads());

    VcfSampleFields sampleFields(referenceSizeInUnits);

    const int shortAlleleSize = genotype.shortAlleleSizeInUnits();
    ReadType shortAlleleSupportType = deterimineSupportType(
        maxUnitsInRead, repeatFindings.countsOfSpanningReads(), repeatFindings.countsOfFlankingReads(),
        shortAlleleSize);

    sampleFields.addAlleleInfo(
        shortAlleleSize, shortAlleleSupportType, genotype.shortAlleleSizeInUnitsCi(),
        readSupportCalculator.getCountOfConsistentSpanningReads(shortAlleleSize),
        readSupportCalculator.getCountOfConsistentFlankingReads(shortAlleleSize),
        readSupportCalculator.getCountOfConsistentRepeatReads(shortAlleleSize));

    if (genotype.numAlleles() == 2)
    {
        const int longAlleleSize = genotype.longAlleleSizeInUnits();
        ReadType longAlleleSupportType = deterimineSupportType(
            maxUnitsInRead, repeatFindings.countsOfSpanningReads(), repeatFindings.countsOfFlankingReads(),
            longAlleleSize);

        sampleFields.addAlleleInfo(
            genotype.longAlleleSizeInUnits(), longAlleleSupportType, genotype.longAlleleSizeInUnitsCi(),
            readSupportCalculator.getCountOfConsistentSpanningReads(longAlleleSize),
            readSupportCalculator.getCountOfConsistentFlankingReads(longAlleleSize),
            readSupportCalculator.getCountOfConsistentRepeatReads(longAlleleSize));
    }

    return sampleFields.encode();
}

void VariantVcfWriter::visit(const RepeatFindings* repeatFindingsPtr)
{
    if (!repeatFindingsPtr->optionalGenotype())
    {
        return;
    }

    const auto& genotype = *(repeatFindingsPtr->optionalGenotype());

    const auto& referenceLocus = variantSpec_.referenceLocus();
    const auto repeatNodeId = variantSpec_.nodes().front();
    const string& repeatUnit = regionSpec_.regionGraph().nodeSeq(repeatNodeId);

    const int referenceSizeInUnits = referenceLocus.length() / repeatUnit.length();

    const string altSymbol = computeAltSymbol(genotype, referenceSizeInUnits);
    const string infoFields = computeInfoFields(variantSpec_, repeatUnit);
    const string sampleFields
        = computeSampleFields(sampleParams_.readLength(), variantSpec_, repeatUnit, *repeatFindingsPtr);

    const int posPreceedingRepeat1based = referenceLocus.start();
    const auto& contigName = reference_.contigInfo().getContigName(referenceLocus.contigIndex());
    const string leftFlankingBase
        = reference_.getSequence(contigName, referenceLocus.start() - 1, referenceLocus.start());

    vector<string> vcfRecordElements
        = { contigName, to_string(posPreceedingRepeat1based), ".",         leftFlankingBase, altSymbol, ".", "PASS",
            infoFields, "GT:SO:REPCN:REPCI:ADSP:ADFL:ADIR",   sampleFields };

    out_ << boost::algorithm::join(vcfRecordElements, "\t") << std::endl;
}

void VariantVcfWriter::visit(const SmallVariantFindings* smallVariantFindingsPtr)
{
    const auto& referenceLocus = variantSpec_.referenceLocus();
    const auto& contigName = reference_.contigInfo().getContigName(referenceLocus.contigIndex());
    string refSequence;
    string altSequence;
    int64_t startPosition = -1;

    if ((variantSpec_.classification().subtype == VariantSubtype::kSwap)
        || (variantSpec_.classification().subtype == VariantSubtype::kSMN))
    {
        assert(variantSpec_.optionalRefNode());
        const auto refNode = *variantSpec_.optionalRefNode();
        const int refNodeIndex = refNode == variantSpec_.nodes().front() ? 0 : 1;
        const int altNodeIndex = refNode == variantSpec_.nodes().front() ? 1 : 0;

        const auto refNodeId = variantSpec_.nodes()[refNodeIndex];
        const auto altNodeId = variantSpec_.nodes()[altNodeIndex];

        refSequence = regionSpec_.regionGraph().nodeSeq(refNodeId);
        altSequence = regionSpec_.regionGraph().nodeSeq(altNodeId);
        // Conversion from 0-based to 1-based coordinates
        startPosition = referenceLocus.start() + 1;
    }
    else if (variantSpec_.classification().subtype == VariantSubtype::kDeletion)
    {
        const string refFlankingBase
            = reference_.getSequence(contigName, referenceLocus.start() - 1, referenceLocus.start());

        const int refNodeId = variantSpec_.nodes().front();
        refSequence = refFlankingBase + regionSpec_.regionGraph().nodeSeq(refNodeId);
        altSequence = refFlankingBase;
        // Conversion from 0-based to 1-based coordinates
        startPosition = referenceLocus.start();
    }
    else if (variantSpec_.classification().subtype == VariantSubtype::kInsertion)
    {
        const string refFlankingBase
            = reference_.getSequence(contigName, referenceLocus.start() - 1, referenceLocus.start());

        const int altNodeId = variantSpec_.nodes().front();
        refSequence = refFlankingBase;
        altSequence = refFlankingBase + regionSpec_.regionGraph().nodeSeq(altNodeId);
        // Conversion from 0-based to 1-based coordinates
        startPosition = referenceLocus.start();
    }
    else
    {
        std::ostringstream encoding;
        encoding << variantSpec_.classification().type << "/" << variantSpec_.classification().subtype;
        throw std::logic_error("Unable to generate VCF record for " + encoding.str());
    }

    const string infoFields = "VARID=" + variantSpec_.id();

    std::vector<string> sampleFields;
    std::vector<string> sampleValues;
    if (variantSpec_.classification().subtype == VariantSubtype::kSMN)
    {
        string genotype;
        switch (smallVariantFindingsPtr->refAllelePresenceStatus())
        {
        case AllelePresenceStatus::kAbsent:
            genotype = "1";
            break;
        case AllelePresenceStatus::kPresent:
            genotype = "0";
            break;
        case AllelePresenceStatus::kUncertain:
            genotype = ".";
            break;
        }
        sampleFields.push_back("GT");
        sampleValues.push_back(genotype);
        sampleFields.push_back("SMN");
        sampleValues.push_back(streamToString(smallVariantFindingsPtr->refAllelePresenceStatus()));
    }
    else
    {
        auto genotype = smallVariantFindingsPtr->optionalGenotype();
        sampleFields.push_back("GT");
        sampleValues.push_back((genotype != boost::none) ? streamToString(*genotype) : ".");
    }

    const string sampleField = boost::algorithm::join(sampleFields, ":");
    const string sampleValue = boost::algorithm::join(sampleValues, ":");

    vector<string> line{
        contigName, streamToString(startPosition), ".", refSequence, altSequence, ".", "PASS", infoFields, sampleField,
        sampleValue
    };
    out_ << boost::algorithm::join(line, "\t") << std::endl;
}

}
