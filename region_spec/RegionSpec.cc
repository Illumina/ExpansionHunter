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

#include "region_spec/RegionSpec.hh"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "thirdparty/json/json.hpp"
#include "thirdparty/spdlog/spdlog.h"

#include "common/common.h"
#include "common/ref_genome.h"

using std::map;
using std::ostream;
using std::string;
using std::to_string;
using std::vector;

using Json = nlohmann::json;

namespace spd = spdlog;

enum class InputRecordType
{
    kRegionWithSingleRepeat,
    kRegionWithMultipleRepeats,
    kUnknown
};

RegionSpec::RegionSpec(
    const string& regionId, const RegionBlueprint& regionBlueprint, AlleleCount expectedAlleleCount,
    const Region& referenceRegion)
    : regionId_(regionId)
    , regionBlueprint_(regionBlueprint)
    , expectedAlleleCount_(expectedAlleleCount)
    , referenceRegion_(referenceRegion)
{
}

bool RegionSpec::operator==(const RegionSpec& other) const
{
    return regionId_ == other.regionId_ && regionBlueprint_ == other.regionBlueprint_
        && offtargetRegions_ == other.offtargetRegions_ && expectedAlleleCount_ == other.expectedAlleleCount_;
}

AlleleCount determineExpectedAlleleCount(Sex sex, const string& chrom)
{
    const bool isFemaleChromY = sex == Sex::kFemale && (chrom == "chrY" || chrom == "Y");
    if (isFemaleChromY)
    {
        return AlleleCount::kZero;
    }

    const bool isSexChrom = chrom == "chrX" || chrom == "X" || chrom == "chrY" || chrom == "Y";
    if (sex == Sex::kMale && isSexChrom)
    {
        return AlleleCount::kOne;
    }

    return AlleleCount::kTwo;
}

// Fill out prefix and suffix sequences
static void loadRegionSequences(
    const RefGenome& reference, const Region& repeatRegion, string& leftFlank, string& repeat, string& rightFlank)
{
    // Reference repeat flanks should be at least as long as reads.
    const int kFlankLen = 1500;

    const int64_t leftFlankStart = repeatRegion.start() - kFlankLen;
    const int64_t leftFlankEnd = repeatRegion.start() - 1;
    const int64_t rightFlankStart = repeatRegion.end() + 1;
    const int64_t rightFlankEnd = repeatRegion.end() + kFlankLen;

    const string leftFlankCoords
        = repeatRegion.chrom() + ":" + to_string(leftFlankStart) + "-" + to_string(leftFlankEnd);
    const string rightFlankCoords
        = repeatRegion.chrom() + ":" + to_string(rightFlankStart) + "-" + to_string(rightFlankEnd);

    reference.ExtractSeq(leftFlankCoords, &leftFlank);
    reference.ExtractSeq(rightFlankCoords, &rightFlank);

    // Output prefix, suffix, repeat, and whole locus.
    const string repeatCoords
        = repeatRegion.chrom() + ":" + to_string(leftFlankEnd + 1) + "-" + to_string(rightFlankStart - 1);
    reference.ExtractSeq(repeatCoords, &repeat);
}

static bool checkIfFieldExists(const Json& record, const string& fieldName)
{
    return record.find(fieldName) != record.end();
}

static void assertFieldExists(const Json& record, const string& fieldName)
{
    if (!checkIfFieldExists(record, fieldName))
    {
        std::stringstream out;
        out << record;
        throw std::logic_error("Field " + fieldName + " must be present in " + out.str());
    }
}

static void assertRecordIsArray(const Json& record)
{
    if (!record.is_array())
    {
        std::stringstream out;
        out << record;
        throw std::logic_error("Expected array but got this instead " + out.str());
    }
}

static InputRecordType guessRecordType(const Json& record)
{
    if (record.find("RepeatId") != record.end())
    {
        return InputRecordType::kRegionWithSingleRepeat;
    }
    else if (record.find("RegionId") != record.end())
    {
        return InputRecordType::kRegionWithMultipleRepeats;
    }
    else
    {
        return InputRecordType::kUnknown;
    }
}

static int getNumberOfRepeats(const string& encoding)
{
    int numBrackets = std::count(encoding.begin(), encoding.end(), '(');
    return numBrackets == 0 ? 1 : numBrackets;
}

static RegionBlueprintComponent::Rarity decodeRarity(const string& encoding)
{
    if (encoding == "common")
    {
        return RegionBlueprintComponent::Rarity::kCommon;
    }
    else if (encoding == "rare")
    {
        return RegionBlueprintComponent::Rarity::kRare;
    }
    else
    {
        throw std::logic_error("Invalid repeat status: " + encoding);
    }
}

static RegionSpec loadSingleRepeatRecord(const Json& record, Sex sampleSex, const RefGenome& reference)
{
    assertFieldExists(record, "RepeatId");
    const string& repeatId = record["RepeatId"];

    assertFieldExists(record, "RepeatUnit");
    const string& repeatUnit = record["RepeatUnit"];

    if (getNumberOfRepeats(repeatUnit) != 1) // TODO: Check that this is valid reference sequence instead
    {
        throw std::runtime_error("Expected repeat unit but got this instead: " + repeatUnit);
    }

    assertFieldExists(record, "ReferenceLocus");
    const string& referenceLocusEncoding = record["ReferenceLocus"];
    const Region referenceLocus(referenceLocusEncoding);

    assertFieldExists(record, "RepeatStatus");
    const string& repeatRarityEncoding = record["RepeatStatus"];

    const auto repeatRarity = decodeRarity(repeatRarityEncoding);

    vector<Region> offtargetRegions;
    if (checkIfFieldExists(record, "OfftargetLoci"))
    {
        assertRecordIsArray(record["OfftargetLoci"]);
        for (const string& locusEncoding : record["OfftargetLoci"])
        {
            offtargetRegions.push_back(Region(locusEncoding));
        }
    }

    string leftFlankSequence, repeatReferenceSequence, rightFlankSequence;
    loadRegionSequences(reference, referenceLocus, leftFlankSequence, repeatReferenceSequence, rightFlankSequence);

    RegionBlueprint blueprint(
        leftFlankSequence, repeatUnit, rightFlankSequence, { repeatId }, { referenceLocus }, { repeatRarity });

    AlleleCount expectedAlleleCount = determineExpectedAlleleCount(sampleSex, referenceLocus.chrom());
    RegionSpec regionSpec(repeatId, blueprint, expectedAlleleCount, referenceLocus);

    regionSpec.setOfftargetRegions(offtargetRegions);

    return regionSpec;
}

static RegionSpec loadMultiRepeatRecord(const Json& record, Sex sampleSex, const RefGenome& reference)
{
    assertFieldExists(record, "RegionId");
    const string& regionId = record["RegionId"];

    assertFieldExists(record, "RegionStructure");
    const string& regionStructure = record["RegionStructure"];
    const std::size_t numRepeats = getNumberOfRepeats(regionStructure);

    assertFieldExists(record, "RepeatIds");
    assertRecordIsArray(record["RepeatIds"]);
    vector<string> repeatIds = record["RepeatIds"];

    assertFieldExists(record, "ReferenceLoci");
    assertRecordIsArray(record["ReferenceLoci"]);
    vector<string> referenceLociEncodings = record["ReferenceLoci"];

    assertFieldExists(record, "RepeatStatuses");
    assertRecordIsArray(record["RepeatStatuses"]);
    vector<string> repeatStatusesEncodings = record["RepeatStatuses"];

    if (repeatIds.size() != numRepeats || referenceLociEncodings.size() != numRepeats
        || repeatStatusesEncodings.size() != numRepeats)
    {
        std::stringstream out;
        out << record;
        throw std::runtime_error("Expected id, locus, and status for each repeat in " + out.str());
    }

    vector<RegionBlueprintComponent::Rarity> repeatRarities;
    for (const auto& repeatRarityEncoding : repeatStatusesEncodings)
    {
        repeatRarities.push_back(decodeRarity(repeatRarityEncoding));
    }

    vector<Region> referenceLoci;
    for (const auto& referenceLocusEncoding : referenceLociEncodings)
    {
        referenceLoci.emplace_back(referenceLocusEncoding);
    }

    const int maxMergeDistance = 150;
    vector<Region> mergedReferenceLoci = merge(referenceLoci, maxMergeDistance);
    if (mergedReferenceLoci.size() != 1)
    {
        std::stringstream out;
        out << record;
        throw std::runtime_error(
            "Expected reference loci to be closer than " + std::to_string(maxMergeDistance)
            + " from one another: " + out.str());
    }

    Region mergedReferenceLocus = mergedReferenceLoci.front();

    string leftFlankSequence, regionReferenceSequence, rightFlankSequence;
    loadRegionSequences(
        reference, mergedReferenceLocus, leftFlankSequence, regionReferenceSequence, rightFlankSequence);

    RegionBlueprint blueprint(
        leftFlankSequence, regionStructure, rightFlankSequence, repeatIds, referenceLoci, repeatRarities);

    AlleleCount expectedAlleleCount = determineExpectedAlleleCount(sampleSex, mergedReferenceLocus.chrom());
    RegionSpec regionSpec(regionId, blueprint, expectedAlleleCount, mergedReferenceLocus);

    // regionSpec.setOfftargetRegions(offtargetRegions);

    return regionSpec;
}

static RegionSpec loadRegionSpecFromJson(const Json& record, Sex sampleSex, const RefGenome& reference)
{
    switch (guessRecordType(record))
    {
    case InputRecordType::kRegionWithSingleRepeat:
        return loadSingleRepeatRecord(record, sampleSex, reference);
    case InputRecordType::kRegionWithMultipleRepeats:
        return loadMultiRepeatRecord(record, sampleSex, reference);
    default:
        std::stringstream out;
        out << record;
        throw std::logic_error("Unknown record type: " + out.str());
    }
}

RegionCatalog loadRegionSpecsFromDisk(const string& specs_path, const RefGenome& reference, Sex sample_sex)
{
    std::ifstream input_stream(specs_path.c_str());

    if (!input_stream.is_open())
    {
        throw std::runtime_error("Failed to open region JSON file " + specs_path);
    }

    Json json_with_region_specs;
    input_stream >> json_with_region_specs;

    if (json_with_region_specs.type() != Json::value_t::array)
    {
        json_with_region_specs = Json::array({ json_with_region_specs });
    }

    RegionCatalog region_specs;

    for (const auto& json_with_region_spec : json_with_region_specs)
    {
        RegionSpec region_spec = loadRegionSpecFromJson(json_with_region_spec, sample_sex, reference);

        // region_spec.setRegionWithInformativeReads(repeat_region.Extend(region_extension_len));

        const string& repeat_id = region_spec.regionId();
        region_specs.emplace(std::make_pair(repeat_id, region_spec));
    }

    return region_specs;
}

ostream& operator<<(ostream& out, const RegionSpec& regionSpec)
{
    for (const auto& component : regionSpec.regionBlueprint())
    {
        out << component;
    }

    return out;
}
