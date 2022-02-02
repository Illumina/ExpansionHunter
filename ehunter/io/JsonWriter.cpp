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

#include "io/JsonWriter.hh"

#include <iomanip>
#include <sstream>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <boost/optional.hpp>

#include "core/Common.hh"
#include "core/ReadSupportCalculator.hh"

namespace ehunter
{

using std::map;
using std::string;
using Json = nlohmann::json;
using boost::optional;
using std::to_string;
using std::vector;

std::ostream& operator<<(std::ostream& out, JsonWriter& jsonWriter)
{
    jsonWriter.write(out);
    return out;
}

JsonWriter::JsonWriter(
    const SampleParameters& sampleParams, const ReferenceContigInfo& contigInfo, const RegionCatalog& regionCatalog,
    const SampleFindings& sampleFindings)
    : sampleParams_(sampleParams)
    , contigInfo_(contigInfo)
    , regionCatalog_(regionCatalog)
    , sampleFindings_(sampleFindings)
{
}

void JsonWriter::write(std::ostream& out)
{
    Json sampleParametersRecord;
    sampleParametersRecord["SampleId"] = sampleParams_.id();
    sampleParametersRecord["Sex"] = streamToString(sampleParams_.sex());

    Json resultsRecord;
    const unsigned locusCount(sampleFindings_.size());
    for (unsigned locusIndex(0); locusIndex < locusCount; ++locusIndex)
    {
        const LocusSpecification& locusSpec = regionCatalog_[locusIndex];
        const LocusFindings& locusFindings = sampleFindings_[locusIndex];
        const std::string& locusId(locusSpec.locusId());

        Json locusRecord;
        locusRecord["LocusId"] = locusId;
        locusRecord["Coverage"] = locusFindings.stats.depth();
        locusRecord["ReadLength"] = locusFindings.stats.meanReadLength();
        locusRecord["FragmentLength"] = locusFindings.stats.medianFragLength();
        locusRecord["AlleleCount"] = static_cast<int>(locusFindings.stats.alleleCount());

        Json variantRecords;
        for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
        {
            const string& variantId = variantIdAndFindings.first;
            const VariantSpecification& variantSpec = locusSpec.getVariantSpecById(variantId);

            VariantJsonWriter variantWriter(contigInfo_, locusSpec, variantSpec);
            variantIdAndFindings.second->accept(&variantWriter);
            variantRecords[variantId] = variantWriter.record();
        }

        if (!variantRecords.empty())
        {
            locusRecord["Variants"] = variantRecords;
        }
        resultsRecord[locusId] = locusRecord;
    }

    Json sampleRecords;
    if (!resultsRecord.empty())
    {
        sampleRecords["LocusResults"] = resultsRecord;
    }
    sampleRecords["SampleParameters"] = sampleParametersRecord;

    out << std::setw(2) << sampleRecords << std::endl;
}

static string encodeGenotype(const RepeatGenotype& genotype)
{
    string encoding = std::to_string(genotype.shortAlleleSizeInUnits());

    if (genotype.numAlleles() == 2)
    {
        encoding = encoding + "/" + std::to_string(genotype.longAlleleSizeInUnits());
    }

    return encoding;
}

void VariantJsonWriter::visit(const RepeatFindings* repeatFindingsPtr)
{
    assert(variantSpec_.classification().type == VariantType::kRepeat);

    const RepeatFindings& repeatFindings = *repeatFindingsPtr;

    record_.clear();
    record_["VariantId"] = variantSpec_.id();
    record_["ReferenceRegion"] = encode(contigInfo_, variantSpec_.referenceLocus());
    record_["VariantType"] = streamToString(variantSpec_.classification().type);
    record_["VariantSubtype"] = streamToString(variantSpec_.classification().subtype);

    const auto repeatNodeId = variantSpec_.nodes().front();
    const auto& repeatUnit = locusSpec_.regionGraph().nodeSeq(repeatNodeId);
    record_["RepeatUnit"] = repeatUnit;

    record_["CountsOfSpanningReads"] = streamToString(repeatFindings.countsOfSpanningReads());
    record_["CountsOfFlankingReads"] = streamToString(repeatFindings.countsOfFlankingReads());
    record_["CountsOfInrepeatReads"] = streamToString(repeatFindings.countsOfInrepeatReads());

    if (repeatFindings.optionalGenotype())
    {
        record_["Genotype"] = encodeGenotype(*repeatFindings.optionalGenotype());
        record_["GenotypeConfidenceInterval"] = streamToString(*repeatFindings.optionalGenotype());
    }

    const auto rfc1Status(repeatFindings.getRFC1Status());
    if (rfc1Status)
    {
        nlohmann::json rfc1Results;
        rfc1Results["Call"] = label(rfc1Status->call);
        rfc1Results["Description"] = rfc1Status->description;
        record_["RFC1MotifAnalysis"] = rfc1Results;
    }
}

void VariantJsonWriter::visit(const SmallVariantFindings* smallVariantFindingsPtr)
{
    const SmallVariantFindings& findings = *smallVariantFindingsPtr;
    record_.clear();
    record_["VariantId"] = variantSpec_.id();
    record_["VariantType"] = streamToString(variantSpec_.classification().type);
    record_["VariantSubtype"] = streamToString(variantSpec_.classification().subtype);
    record_["ReferenceRegion"] = encode(contigInfo_, variantSpec_.referenceLocus());
    record_["CountOfRefReads"] = findings.numRefReads();
    record_["CountOfAltReads"] = findings.numAltReads();
    record_["StatusOfRefAllele"] = streamToString(findings.refAllelePresenceStatus().status);
    record_["LogLikelihoodRefAllelePresent"] = streamToString(findings.refAllelePresenceStatus().logLikelihoodRatio);
    record_["StatusOfAltAllele"] = streamToString(findings.altAllelePresenceStatus().status);
    record_["LogLikelihoodAltAllelePresent"] = streamToString(findings.altAllelePresenceStatus().logLikelihoodRatio);
    if (findings.optionalGenotype())
    {
        record_["Genotype"] = streamToString(*findings.optionalGenotype());
    }
}

}
