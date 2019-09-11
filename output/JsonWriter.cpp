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

#include "output/JsonWriter.hh"

#include <iomanip>
#include <sstream>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <boost/optional.hpp>

#include "common/Common.hh"
#include "stats/ReadSupportCalculator.hh"

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
    for (const auto& locusIdAndFindings : sampleFindings_)
    {
        const string& locusId = locusIdAndFindings.first;
        const LocusSpecification& locusSpec = regionCatalog_.at(locusId);
        const LocusFindings& locusFindings = locusIdAndFindings.second;

        Json locusRecord;
        locusRecord["LocusId"] = locusId;

        if (locusFindings.optionalStats)
        {
            locusRecord["AlleleCount"] = static_cast<int>(locusFindings.optionalStats->alleleCount());
            locusRecord["Coverage"] = locusFindings.optionalStats->depth();
            locusRecord["ReadLength"] = locusFindings.optionalStats->meanReadLength();
        }

        Json variantRecords;
        for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
        {
            const string& variantId = variantIdAndFindings.first;
            const VariantSpecification& variantSpec = locusSpec.getVariantSpecById(variantId);

            VariantJsonWriter variantWriter(contigInfo_, locusSpec, variantSpec);
            variantIdAndFindings.second->accept(variantWriter);
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

void VariantJsonWriter::visit(StrFindings& strFindings)
{
    assert(variantSpec_.classification().type == VariantType::kRepeat);

    record_.clear();
    record_["VariantId"] = variantSpec_.id();
    record_["ReferenceRegion"] = encode(contigInfo_, variantSpec_.referenceLocus());
    record_["VariantType"] = streamToString(variantSpec_.classification().type);
    record_["VariantSubtype"] = streamToString(variantSpec_.classification().subtype);

    const auto repeatNodeId = variantSpec_.nodes().front();
    const auto& repeatUnit = locusSpec_.regionGraph().nodeSeq(repeatNodeId);
    record_["RepeatUnit"] = repeatUnit;

    record_["CountsOfSpanningReads"] = streamToString(strFindings.countsOfSpanningReads());
    record_["CountsOfFlankingReads"] = streamToString(strFindings.countsOfFlankingReads());
    record_["CountsOfInrepeatReads"] = streamToString(strFindings.countsOfInrepeatReads());

    if (strFindings.optionalGenotype())
    {
        record_["Genotype"] = encodeGenotype(*strFindings.optionalGenotype());
        record_["GenotypeConfidenceInterval"] = streamToString(*strFindings.optionalGenotype());
    }
}

/*
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
} */

}
