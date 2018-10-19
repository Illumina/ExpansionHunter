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

// clang-format off
const string kVcfHeaderBase =
    "##fileformat=VCFv4.1\n"
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n"
    "##INFO=<ID=REF,Number=1,Type=Integer,Description=\"Reference copy number\">\n"
    "##INFO=<ID=RL,Number=1,Type=Integer,Description=\"Reference length in bp\">\n"
    "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat unit in the reference orientation\">\n"
    "##INFO=<ID=REPID,Number=1,Type=String,Description=\"Repeat identifier from the repeat specification file\">\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "##FORMAT=<ID=SO,Number=1,Type=String,Description=\"Type of reads that support the allele; can be SPANNING, "
        "FLANKING, or INREPEAT meaning that the reads span, flank, or are fully contained in the repeat\">\n"
    "##FORMAT=<ID=REPCN,Number=1,Type=String,Description=\"Number of repeat units spanned by the allele\">\n"
    "##FORMAT=<ID=REPCI,Number=1,Type=String,Description=\"Confidence interval for REPCN\">\n"
    "##FORMAT=<ID=ADFL,Number=1,Type=String,Description=\"Number of flanking reads consistent with the allele\">\n"
    "##FORMAT=<ID=ADSP,Number=1,Type=String,Description=\"Number of spanning reads consistent with the allele\">\n"
    "##FORMAT=<ID=ADIR,Number=1,Type=String,Description=\"Number of in-repeat reads consistent with the allele\">\n";
// clang-format on

void writeRepeatSizeHeaderLine(int repeatSize, ostream& out)
{
    out << "##ALT=<ID=STR" << repeatSize << ",Description=\"Allele comprised of " << repeatSize << " repeat units\">\n";
}

void writeBodyHeader(const string& sampleName, ostream& out)
{
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sampleName << "\n";
}

std::ostream& operator<<(std::ostream& out, VcfWriter& vcfWriter)
{
    vcfWriter.writeHeader(out);
    vcfWriter.writeBody(out);
    return out;
}

VcfWriter::VcfWriter(
    const string& sampleName, int readLength, const RegionCatalog& regionSpecs, const SampleFindings& sampleFindings)
    : sampleName_(sampleName)
    , readLength_(readLength)
    , regionSpecs_(regionSpecs)
    , sampleFindings_(sampleFindings)
{
}

void VcfWriter::writeHeader(ostream& out)
{
    out << kVcfHeaderBase;

    const set<int> altRepeatSizes = computeAltRepeatSizes(regionSpecs_, sampleFindings_);

    for (int repeatSize : altRepeatSizes)
    {
        writeRepeatSizeHeaderLine(repeatSize, out);
    }

    writeBodyHeader(sampleName_, out);
}

void VcfWriter::writeBody(ostream& out)
{
    for (const auto& regionIdAndRegionFindings : sampleFindings_)
    {
        const string& regionId = regionIdAndRegionFindings.first;
        const RegionFindings& regionFindings = regionIdAndRegionFindings.second;

        writeRegionFindings(regionId, regionFindings, out);
    }
}

void VcfWriter::writeRegionFindings(const string& regionId, const RegionFindings& regionFindings, ostream& out)
{
    const RegionSpec& regionSpec = regionSpecs_.at(regionId);

    for (const auto& blueprintComponent : regionSpec.regionBlueprint())
    {
        if (blueprintComponent.type() == RegionBlueprintComponent::Type::kRepeat)
        {
            const auto& repeatIdAndRepeatFindingsIter = regionFindings.find(blueprintComponent.id());
            assert(repeatIdAndRepeatFindingsIter != regionFindings.end());
            const RepeatFindings& repeatFindings = repeatIdAndRepeatFindingsIter->second;

            writeRepeatFindings(blueprintComponent, repeatFindings, out);
        }
    }
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

static string computeInfoFields(const RegionBlueprintComponent& repeatBlueprint)
{
    const string& repeatUnit = repeatBlueprint.sequence();
    const auto referenceRegion = repeatBlueprint.referenceRegion();
    assert(referenceRegion);
    const int referenceSizeInBp = referenceRegion->length();
    const int referenceSizeInUnits = referenceSizeInBp / repeatUnit.length();

    vector<string> fields;
    fields.push_back("END=" + std::to_string(referenceRegion->end()));
    fields.push_back("REF=" + std::to_string(referenceSizeInUnits));
    fields.push_back("RL=" + std::to_string(referenceSizeInUnits));
    fields.push_back("RU=" + repeatUnit);
    fields.push_back("REPID=" + repeatBlueprint.id());

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
    int readLength, const RegionBlueprintComponent& repeatBlueprint, const RepeatFindings& repeatFindings)
{
    const string& repeatUnit = repeatBlueprint.sequence();
    const auto referenceRegion = repeatBlueprint.referenceRegion();
    assert(referenceRegion);
    const int referenceSizeInBp = referenceRegion->length();
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

void VcfWriter::writeRepeatFindings(
    const RegionBlueprintComponent& repeatBlueprint, const RepeatFindings& repeatFindings, ostream& out)
{
    if (!repeatFindings.optionalGenotype())
    {
        return;
    }

    const auto referenceRegion = repeatBlueprint.referenceRegion();
    assert(referenceRegion);

    const char leftFlankingBase = '?';
    const string& repeatUnit = repeatBlueprint.sequence();
    const int referenceSizeInUnits = referenceRegion->length() / repeatUnit.length();

    const string altSymbol = computeAltSymbol(*repeatFindings.optionalGenotype(), referenceSizeInUnits);
    const string infoFields = computeInfoFields(repeatBlueprint);
    const string sampleFields = computeSampleFields(readLength_, repeatBlueprint, repeatFindings);

    out << referenceRegion->chrom() << "\t" << referenceRegion->start() - 1 << "\t.\t" << leftFlankingBase << "\t"
        << altSymbol << "\t.\tPASS\t" << infoFields << "\tGT:SO:REPCN:REPCI:ADSP:ADFL:ADIR\t" << sampleFields
        << std::endl;
}
