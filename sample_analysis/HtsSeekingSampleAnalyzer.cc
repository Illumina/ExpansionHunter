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

#include "sample_analysis/HtsSeekingSampleAnalyzer.hh"

#include <cassert>
#include <memory>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

#include "reads/read_pairs.h"
#include "region_analysis/RegionAnalyzer.hh"
#include "sample_analysis/HtsFileSeeker.hh"
#include "sample_analysis/MateExtractor.hh"

using boost::optional;
using htshelpers::HtsFileSeeker;
using reads::LinearAlignmentStats;
using reads::Read;
using reads::ReadPairs;
using std::ostream;
using std::string;
using std::unordered_map;
using std::vector;

using AlignmentStatsCatalog = unordered_map<string, LinearAlignmentStats>;

static ReadPairs collectRelevantReads(
    int readLength, const RegionSpec& regionSpec, AlignmentStatsCatalog& alignmentStatsCatalog,
    HtsFileSeeker& fileHopper)
{
    ReadPairs readPairs;

    const auto extractionRegion = regionSpec.referenceRegion().Extend(2 * readLength);
    fileHopper.setRegion(extractionRegion);
    while (fileHopper.trySeekingToNextPrimaryAlignment())
    {
        LinearAlignmentStats alignmentStats;
        Read read = fileHopper.decodeRead(alignmentStats);
        readPairs.Add(read);
        alignmentStatsCatalog.emplace(std::make_pair(read.readId(), alignmentStats));
    }

    return readPairs;
}

bool checkIfMatesWereMappedNearby(const LinearAlignmentStats& alignmentStats)
{
    if ((alignmentStats.chrom_id == alignmentStats.mate_chrom_id)
        && (std::abs(alignmentStats.pos - alignmentStats.mate_pos) < 1000))
    {
        return true;
    }
    return false;
}

void recoverMates(const string& htsFilePath, const AlignmentStatsCatalog& alignmentStatsCatalog, ReadPairs& readPairs)
{
    htshelpers::MateExtractor mateExtractor(htsFilePath);

    for (auto& fragmentIdAndReadPair : readPairs)
    {
        reads::ReadPair& readPair = fragmentIdAndReadPair.second;

        if (readPair.first_mate.isSet() && readPair.second_mate.isSet())
        {
            continue;
        }

        assert(readPair.first_mate.isSet() || readPair.second_mate.isSet());

        const Read& read = readPair.first_mate.isSet() ? readPair.first_mate : readPair.second_mate;

        const auto alignmentStatsIterator = alignmentStatsCatalog.find(read.readId());
        if (alignmentStatsIterator == alignmentStatsCatalog.end())
        {
            std::logic_error("Cannon recover mate of uncatalogued read");
        }
        const LinearAlignmentStats& alignmentStats = alignmentStatsIterator->second;

        if (!checkIfMatesWereMappedNearby(alignmentStats))
        {
            Read mate = mateExtractor.extractMate(read, alignmentStats);
            if (!readPair.first_mate.isSet())
            {
                readPair.first_mate = mate;
            }
            else
            {
                readPair.second_mate = mate;
            }
        }
    }
}

static RegionFindings analyzeRegion(
    const ReadPairs& readPairs, const RegionSpec& regionSpec, double haplotypeDepth, int readLength,
    const string& alignerName, ostream& alignmentStream)
{
    alignmentStream << regionSpec.regionId() << ":" << std::endl;
    GraphAlignmentHeuristicsParameters alignmentParams;
    RegionAnalyzer regionAnalyzer(
        regionSpec, haplotypeDepth, readLength, alignmentStream, alignerName, alignmentParams);

    for (const auto fragmentIdAndReads : readPairs)
    {
        const auto& readPair = fragmentIdAndReads.second;
        if (readPair.first_mate.isSet() && readPair.second_mate.isSet())
        {
            regionAnalyzer.processMates(readPair.first_mate, readPair.second_mate);
        }
    }

    return regionAnalyzer.genotype();
}

SampleFindings htsSeekingSampleAnalysis(
    const string& htsFilePath, double haplotypeDepth, int readLength, const RegionCatalog& regionCatalog,
    const string& alignerName, ostream& alignmentStream)
{
    HtsFileSeeker HtsFileSeeker(htsFilePath);

    SampleFindings sampleFindings;
    for (const auto& regionIdAndRegionSpec : regionCatalog)
    {
        const string& regionId = regionIdAndRegionSpec.first;
        const RegionSpec& regionSpec = regionIdAndRegionSpec.second;
        AlignmentStatsCatalog alignmentStatsCatalog;
        ReadPairs readPairs = collectRelevantReads(readLength, regionSpec, alignmentStatsCatalog, HtsFileSeeker);
        recoverMates(htsFilePath, alignmentStatsCatalog, readPairs);

        auto regionFindings
            = analyzeRegion(readPairs, regionSpec, haplotypeDepth, readLength, alignerName, alignmentStream);
        sampleFindings.emplace(std::make_pair(regionId, regionFindings));
    }

    return sampleFindings;
}
