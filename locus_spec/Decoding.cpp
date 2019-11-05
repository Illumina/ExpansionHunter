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

#include "locus_spec/Decoding.hh"

#include "common/WorkflowContext.hh"

namespace ehunter
{

using std::string;
using std::to_string;
using std::vector;

static vector<GenomicRegion> decodeRegions(const Reference& reference, const vector<string>& encodings)
{
    const ReferenceContigInfo& contigInfo = reference.contigInfo();

    vector<GenomicRegion> regions;
    regions.reserve(encodings.size());
    for (const auto& encoding : encodings)
    {
        regions.push_back(decode(contigInfo, encoding));
    }
    return regions;
}

GraphLocusSpec decode(const Reference& reference, const GraphLocusEncoding& locusEncoding)
{

    /*

    vector<GenomicRegion> targetReadExtractionRegions;
    for (const GenomicRegion& region : userDescription.targetRegions)
    {
        targetReadExtractionRegions.push_back(region.extend(kExtensionLength));
    }
    if (targetReadExtractionRegions.empty())
    {

        targetReadExtractionRegions.push_back(userDescription.locusLocation.extend(kExtensionLength));
    } */

    /*



    NodeLocations referenceRegionsOfGraphNodes
        = associateNodesWithReferenceRegions(blueprint, locusGraph, completeReferenceRegions);


    GraphLocusReferenceRegions referenceRegions;
    referenceRegions.offtargetRegionsWithReads = userDescription.offtargetRegions;
    referenceRegions.regionsWithReads = targetReadExtractionRegions;
    referenceRegions.statsRegions = computeStatsRegions(userDescription.locusLocation, kExtensionLength);

    ReferenceGraph referenceGraph(locusGraph, referenceRegionsOfGraphNodes);

    GraphLocusSpec locusSpec(locusId, copyNumber, referenceRegions, referenceGraph, parameters);

    int variantIndex = 0;
    for (const auto& feature : blueprint)
    {
        if (doesFeatureDefineVariant(feature.type))
        {
            VariantDescriptionFromUser variant = userDescription.variantDescriptionFromUsers.at(variantIndex);
            const GenomicRegion& referenceRegion = variant.variantLocation;

            VariantTypeFromUser variantTypeFromDescription = variant.variantType;
            const string& variantId = variant.variantId;
            VariantType variantType = determineVariantType(feature.type);
            VariantSubtype variantSubtype
                = determineVariantSubtype(feature.type, variantTypeFromDescription, referenceRegion);

            optional<NodeId> optionalReferenceNode = determineReferenceNode(feature, reference, referenceRegion);

            VariantClassification classification(variantType, variantSubtype);

            locusSpec.addVariantSpecification(
                variantId, classification, referenceRegion, feature.nodeIds, optionalReferenceNode);

            ++variantIndex;
        }
    } */

    return locusSpec;
}

}