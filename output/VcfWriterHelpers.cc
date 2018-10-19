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

#include "output/VcfWriterHelpers.hh"

#include <boost/algorithm/string/join.hpp>

#include <sstream>
#include <vector>

using std::deque;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::stringstream;
using std::to_string;
using std::vector;

template <typename T> string encodeSampleFields(const deque<T>& fieldRecords)
{
    vector<string> fieldRecordEncodings;

    for (const auto& record : fieldRecords)
    {
        stringstream stringStream;
        stringStream << record;
        fieldRecordEncodings.emplace_back(stringStream.str());
    }

    return boost::algorithm::join(fieldRecordEncodings, "/");
}

string VcfSampleFields::encode() const
{
    vector<string> encoding;
    encoding.push_back(encodeSampleFields(genotype_));
    encoding.push_back(encodeSampleFields(sources_));
    encoding.push_back(encodeSampleFields(alleleSizes_));
    encoding.push_back(encodeSampleFields(confidenceIntervals_));
    encoding.push_back(encodeSampleFields(spanningReadCounts_));
    encoding.push_back(encodeSampleFields(flankingReadCounts_));
    encoding.push_back(encodeSampleFields(repeatReadCounts_));

    return boost::algorithm::join(encoding, ":");
}

void VcfSampleFields::addAlleleInfo(
    int alleleSize, ReadType source, NumericInterval confidenceInterval, int spanningReadCount, int flankingReadCount,
    int repeatReadCount)
{
    if (alleleSize == referenceSize_)
    {
        addRefAlleleInfo(alleleSize, source, confidenceInterval, spanningReadCount, flankingReadCount, repeatReadCount);
    }
    else
    {
        addAltAlleleInfo(alleleSize, source, confidenceInterval, spanningReadCount, flankingReadCount, repeatReadCount);
    }
}

void VcfSampleFields::addRefAlleleInfo(
    int alleleSize, ReadType source, NumericInterval confidenceInterval, int spanningReadCount, int flankingReadCount,
    int repeatReadCount)
{
    genotype_.push_front(0);
    alleleSizes_.push_front(alleleSize);
    sources_.push_front(source);
    confidenceIntervals_.push_front(confidenceInterval);
    spanningReadCounts_.push_front(spanningReadCount);
    flankingReadCounts_.push_front(flankingReadCount);
    repeatReadCounts_.push_front(repeatReadCount);
}

void VcfSampleFields::addAltAlleleInfo(
    int alleleSize, ReadType source, NumericInterval confidenceInterval, int spanningReadCount, int flankingReadCount,
    int repeatReadCount)
{
    const int previousAlleleSize = genotype_.empty() ? -1 : genotype_.back();
    if (alleleSize <= previousAlleleSize)
    {
        throw std::logic_error(
            "Allele of size " + to_string(alleleSize) + " cannot follow allele of size "
            + to_string(previousAlleleSize));
    }

    int haplotypeNum = 1;
    if (!genotype_.empty())
    {
        const int previousAlleleSize = alleleSizes_.back();
        haplotypeNum = previousAlleleSize == alleleSize ? genotype_.back() : genotype_.back() + 1;
    }

    genotype_.push_back(haplotypeNum);
    alleleSizes_.push_back(alleleSize);
    sources_.push_back(source);
    confidenceIntervals_.push_back(confidenceInterval);
    spanningReadCounts_.push_back(spanningReadCount);
    flankingReadCounts_.push_back(flankingReadCount);
    repeatReadCounts_.push_back(repeatReadCount);
}

set<int> computeAltRepeatSizes(const RegionCatalog& regionSpecs, const SampleFindings& sampleFindings)
{
    set<int> altRepeatSizes;

    for (const auto& regionIdAndRegionFindings : sampleFindings)
    {
        const string& regionId = regionIdAndRegionFindings.first;
        const RegionFindings& regionFindings = regionIdAndRegionFindings.second;

        const RegionSpec& regionSpec = regionSpecs.at(regionId);

        for (const auto& blueprintComponent : regionSpec.regionBlueprint())
        {
            if (blueprintComponent.type() != RegionBlueprintComponent::Type::kRepeat)
            {
                continue;
            }

            const auto& repeatIdAndRepeatFindingsIter = regionFindings.find(blueprintComponent.id());
            assert(repeatIdAndRepeatFindingsIter != regionFindings.end());
            const RepeatFindings& repeatFindings = repeatIdAndRepeatFindingsIter->second;

            const string& repeatUnit = blueprintComponent.sequence();
            const auto referenceRegion = blueprintComponent.referenceRegion();
            assert(referenceRegion);
            const int referenceSize = referenceRegion->length() / repeatUnit.length();

            if (!repeatFindings.optionalGenotype())
            {
                continue;
            }

            const RepeatGenotype& genotype = repeatFindings.optionalGenotype().get();

            if (genotype.shortAlleleSizeInUnits() != referenceSize)
            {
                altRepeatSizes.insert(genotype.shortAlleleSizeInUnits());
            }

            if (genotype.longAlleleSizeInUnits() != referenceSize)
            {
                altRepeatSizes.insert(genotype.longAlleleSizeInUnits());
            }
        }
    }

    return altRepeatSizes;
}
