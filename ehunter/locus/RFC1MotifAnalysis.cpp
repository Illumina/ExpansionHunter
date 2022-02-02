//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Chris Saunders <csaunders@illumina.com>
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

#include "RFC1MotifAnalysis.hh"

#include <algorithm>
#include <numeric>

#include "boost/algorithm/string.hpp"
#include "boost/optional.hpp"
#include "boost/range/adaptor/map.hpp"
#include "boost/range/numeric.hpp"

#include "graphalign/GraphAlignment.hh"
#include "locus/AlignmentBuffer.hh"
#include "locus/RFC1MotifAnalysisUtil.hh"
#include "locus/RFC1Status.hh"

namespace ehunter
{

namespace
{

/// \brief Get the number of tandem repeat motifs spanned by the given alignment
///
/// This assumes that group 1 is the repeat motif -- the hard coded group number works for RFC1 but isn't a general
/// repeat motif counter for other loci.
///
unsigned getRepeatMotifSpan(const graphtools::GraphAlignment& readAlign)
{
    unsigned count(0);
    for (size_t index(0); index < readAlign.size(); ++index)
    {
        if (readAlign.path().getNodeIdByIndex(index) == 1)
        {
            count++;
        }
    }
    return count;
}

/// \brief Analyze the read alignment to get the total aligned length to each graph node
///
/// Any node which the read does not align to will have a length of zero
///
std::vector<unsigned> getGraphNodeAlignmentLengths(const graphtools::GraphAlignment& readAlign)
{
    std::vector<unsigned> nodeLengths;
    for (size_t index(0); index < readAlign.size(); ++index)
    {
        const auto node(readAlign.path().getNodeIdByIndex(index));
        if (node + 1 > nodeLengths.size())
        {
            nodeLengths.resize(node + 1, 0);
        }
        nodeLengths[node] += readAlign[index].queryLength();
    }

    return nodeLengths;
}

/// \brief Extract repeat motif genotype value from repeat findings
///
std::vector<unsigned> getAlleleRepeatMotifCounts(const RepeatFindings& repeatFindings)
{
    std::vector<unsigned> counts;
    assert(repeatFindings.optionalGenotype());
    const RepeatGenotype& genotype(*repeatFindings.optionalGenotype());

    counts.push_back(genotype.shortAlleleSizeInUnits());
    if (genotype.numAlleles() == 2)
    {
        counts.push_back(genotype.longAlleleSizeInUnits());
    }
    return counts;
}

/// \brief Extract basecall quality information from EH read string
///
/// In EH, the read string encodes a binary low and high quality state using case. This routine extracts a binary
/// quality string from the read.
///
std::vector<uint8_t> getBinaryQuals(const std::string& read)
{
    std::vector<uint8_t> binaryQuals(read.size());
    std::transform(read.begin(), read.end(), binaryQuals.begin(), [](unsigned char c) { return std::isupper(c); });
    return binaryQuals;
}

using MotifsWithQualWeight_t = std::vector<std::pair<std::string, double>>;

/// \brief Given a single EH alignment record, return all repeat motifs with quality information
///
/// All motifs are found and reported in the read's alignment orientation
///
/// \param[in] binaryQuals Quality vector for the read, reduced to 2 {low, high} quality states
///
/// \param[in] usableBaseRange 0-indexed, closed interval representing the usable portion of the read in read alignment
/// coordinates, or none. If none, the entire read is considered usable.
///
/// \param[in] expectedMotifSize Size of the repeat motif at the target locus
///
/// \return A vector of motif data. Each vector element is a 2-tuple containing the motif string and its quality weight.
/// The quality weight is the fraction of high-quality bases in the motif.
///
MotifsWithQualWeight_t getAllMotifsWithQualWeight(
    const locus::AlignmentBufferData& alignmentData, const std::vector<uint8_t>& binaryQuals,
    boost::optional<std::pair<unsigned, unsigned>> usableBaseRange, const unsigned expectedMotifSize)
{
    // A motif most be found at least this far from the edge of a read for it to be counted
    static const unsigned minDistFromReadEdge(1);

    const unsigned readLength(binaryQuals.size());

    // Get the first and last positions in read coordinates from which a repeat motif can be extracted
    // (zero-indexed, closed)
    const auto repeatExtractionRange
        = std::make_pair(minDistFromReadEdge, readLength - (expectedMotifSize + minDistFromReadEdge));

    auto isValidMotifCandidate = [&](const unsigned readPos) -> bool
    {
        // Exclude motifs in low-quality regions of the read:
        if (usableBaseRange)
        {
            if ((readPos < usableBaseRange->first) or ((readPos + expectedMotifSize - 1) > usableBaseRange->second))
            {
                return false;
            }
        }

        // Skip repeats occurring at the start or end of the read
        if ((readPos < repeatExtractionRange.first) or (readPos > repeatExtractionRange.second))
        {
            return false;
        }

        return true;
    };

    const auto read(boost::to_upper_copy<std::string>(alignmentData.read));

    // Get the first and last usable positions in read coordinates of the repeat tract (zero-indexed, closed)
    //
    // Note that the first and last repeat units in the tract are not used (presumably to reduce motif noise?)
    auto repeatTractRange(std::make_pair(0, readLength - 1));

    const auto nodeAlignmentLengths = getGraphNodeAlignmentLengths(alignmentData.readAlignment);
    if ((nodeAlignmentLengths.size() > 0) and (nodeAlignmentLengths[0] != 0))
    {
        repeatTractRange.first += nodeAlignmentLengths[0] + expectedMotifSize;
    }
    if ((nodeAlignmentLengths.size() > 2) and (nodeAlignmentLengths[2] != 0))
    {
        repeatTractRange.second -= nodeAlignmentLengths[2] + expectedMotifSize;
    }

    std::vector<std::pair<std::string, double>> motifData;
    for (unsigned readPos(repeatTractRange.first); readPos < (repeatTractRange.second + 1);
         readPos += expectedMotifSize)
    {
        if (isValidMotifCandidate(readPos))
        {
            const auto motifStr(read.substr(readPos, expectedMotifSize));
            const auto qiter(binaryQuals.begin() + readPos);
            const double motifQualWeight = mean(qiter, qiter + expectedMotifSize);
            motifData.emplace_back(getMinRotation(motifStr), motifQualWeight);
        }
    }

    return motifData;
}

using MotifCountMap_t = std::map<std::string, unsigned>;

/// \brief Iterate over reads from the RFC1 locus and count high quality repeat motif observations
///
/// \param[in] expectedMotifSize Size of the repeat motif at the target locus
///
/// \param[in] minRepeatMotifSpan The minimum number of repeat units which must be spanned by a read alignment for the
/// read to be used as evidence
///
/// \return A map from high-quality motif strings to the number of times the motif has been observed in the target STR
/// region
///
MotifCountMap_t getHighQMotifMap(
    const locus::AlignmentBuffer& alignmentBuffer, const unsigned expectedMotifSize, const unsigned minRepatMotifSpan)
{
    // Only include motifs with a motif quality weight at least this high, where the motif quality weight is the
    // fraction of motif bases that are quantized as the "high quality" state.
    static const double minMotifQualWeight(1.0);

    MotifCountMap_t highQMotifMap;
    for (const auto& alignmentData : alignmentBuffer.getBuffer())
    {
        // Only look at reads which align to the repeat unit at least minRepeatMotifSpan times:
        if (getRepeatMotifSpan(alignmentData.readAlignment) < minRepatMotifSpan)
        {
            continue;
        }

        const auto binaryQuals = getBinaryQuals(alignmentData.read);
        const auto usableBaseRange = findUsableReadBaseRange(binaryQuals, alignmentData.isReversed);
        if (not usableBaseRange)
        {
            continue;
        }

        const auto motifData
            = getAllMotifsWithQualWeight(alignmentData, binaryQuals, usableBaseRange, expectedMotifSize);

        for (const auto& elem : motifData)
        {
            if (elem.second >= minMotifQualWeight)
            {
                highQMotifMap[elem.first] += 1;
            }
        }
    }

    return highQMotifMap;
}

/// \brief Observation data for a given motif
///
struct MotifObservations
{
    /// Total motif observation count
    unsigned count = 0;

    /// Total quality-weighted motif observation count
    double weightedCount = 0.;

    /// Fraction of total quality-weighted motif observation count over all motifs
    double weightedFrac = 0.;
};

using MotifObservationMap_t = std::map<std::string, MotifObservations>;

/// \brief Summary data on all high-quality repeat motif observations, as well as the fraction of pathogenic motifs per
/// read.
///
struct MotifAndPurityData
{
    /// A map from high quality repeat motifs to motif observation data
    MotifObservationMap_t motifMap;

    /// For each read, contains the fraction of pathogenic motifs compared to all other high-quality repeat motifs in
    /// the read
    std::vector<double> pathogenicMotifFractionPerRead;
};

using HighQMotif_t = std::set<std::string>;

/// \brief Find set of high quality motifs at the RFC1 locus
///
/// \param[in] expectedMotifSize: Size of the repeat motif at the target locus
///
/// \param[in] minRepeatMotifSpan The minimum number of repeat units which must be spanned by a read alignment for the
/// read to be used as evidence
///
HighQMotif_t gethighQMotifs(
    const locus::AlignmentBuffer& alignmentBuffer, const unsigned expectedMotifSize, const unsigned minRepeatMotifSpan)
{
    // A motif must have at least this many high-quality observations before it is included in the highQ motif set
    static const unsigned minHighQMotifObservations(2);

    const auto highQMotifMap = getHighQMotifMap(alignmentBuffer, expectedMotifSize, minRepeatMotifSpan);

    HighQMotif_t highQMotifs;
    for (const auto& elem : highQMotifMap)
    {
        const auto& highQMotif(elem.first);
        const auto count(elem.second);
        if ((highQMotif.size() == expectedMotifSize) and (count >= minHighQMotifObservations))
        {
            highQMotifs.insert(highQMotif);
        }
    }

    return highQMotifs;
}

/// \brief Return the total pathogenic motif observations
///
/// \param[in] readMotifCount Count of motif observations
///
/// \param[in] pathogenicMotifs Container of pathogenic motif strings
///
unsigned getPathogenicMotifTotal(
    const std::map<std::string, unsigned> readMotifCount, const std::vector<std::string>& pathogenicMotifs)
{
    unsigned total(0);
    for (const auto& pathogenicMotif : pathogenicMotifs)
    {
        const auto citer(readMotifCount.find(pathogenicMotif));
        if (citer != readMotifCount.end())
        {
            total += citer->second;
        }
    }
    return total;
}

/// \brief Get MotifAndPurityData from RFC1 locus reads
///
/// \param[in] expectedMotifSize: Size of the repeat motif at the target locus
///
/// \param[in] minRepeatMotifSpan The minimum number of repeat units which must be spanned by a read alignment for the
/// read to be used as evidence
///
/// \param[in] pathogenicMotifs Container of pathogenic motif strings
///
MotifAndPurityData getMotifAndPurityData(
    const locus::AlignmentBuffer& alignmentBuffer, const unsigned expectedMotifSize, const unsigned minRepatMotifSpan,
    const std::vector<std::string>& pathogenicMotifs)
{
    // For a read to be counted in the pathogen_purities list, at least this many repeat motifs must be processed from
    // the read alignment
    static const unsigned minPurityMotifCountsPerRead(5);

    const auto highQMotifs(gethighQMotifs(alignmentBuffer, expectedMotifSize, minRepatMotifSpan));

    MotifAndPurityData mpData;

    for (const auto& alignmentData : alignmentBuffer.getBuffer())
    {
        // Only look at reads which align to the repeat unit at least minRepeatMotifSpan times:
        if (getRepeatMotifSpan(alignmentData.readAlignment) < minRepatMotifSpan)
        {
            continue;
        }

        const auto binaryQuals = getBinaryQuals(alignmentData.read);
        const auto motifData = getAllMotifsWithQualWeight(alignmentData, binaryQuals, {}, expectedMotifSize);

        std::map<std::string, unsigned> readMotifCount;
        for (const auto& elem : motifData)
        {
            const auto& motifSeq(elem.first);
            const auto motifQualWeight(elem.second);
            if (highQMotifs.find(motifSeq) != highQMotifs.end())
            {
                mpData.motifMap[motifSeq].count += 1;
                mpData.motifMap[motifSeq].weightedCount += motifQualWeight;
                readMotifCount[motifSeq] += 1;
            }
        }

        const unsigned readMotifTotal = boost::accumulate(readMotifCount | boost::adaptors::map_values, 0u);
        if (readMotifTotal >= minPurityMotifCountsPerRead)
        {
            const unsigned pathogenicMotifTotal(getPathogenicMotifTotal(readMotifCount, pathogenicMotifs));
            if (pathogenicMotifTotal > 0)
            {
                mpData.pathogenicMotifFractionPerRead.push_back(safeFrac(pathogenicMotifTotal, readMotifTotal));
            }
        }
    }

    double totalWeightedCount(0);
    for (auto& v : mpData.motifMap | boost::adaptors::map_values)
    {
        totalWeightedCount += v.weightedCount;
    }
    for (auto& v : mpData.motifMap | boost::adaptors::map_values)
    {
        v.weightedFrac = v.weightedCount / totalWeightedCount;
    }

    return mpData;
}

/// \brief Count number of spanning reads for RFC1 locus.
///
/// To be counted, spanning reads must overlap with both left and right flanks by a sufficient amount and also meet
/// quality criteria.
///
/// \return Spanning read count
///
unsigned countSpanningReads(const locus::AlignmentBuffer& alignmentBuffer)
{
    // To count as spanning reads, each flank alignment should be at least this long
    static const int minFlankLength(10);

    // To count as spanning reads, each flank alignment should have at least this fraction of high quality bases
    //
    // This is tested over twice the minFlankLength
    static const double minFlankHighQBaseFraction(0.7);

    unsigned numSpanningReads(0);
    for (const auto& alignmentData : alignmentBuffer.getBuffer())
    {
        const auto nodeAlignmentLengths = getGraphNodeAlignmentLengths(alignmentData.readAlignment);

        // Require that the read aligns to both flanks and the repeat:
        if ((nodeAlignmentLengths.size() < 3)
            or (std::any_of(
                nodeAlignmentLengths.begin(), nodeAlignmentLengths.end(), [](unsigned i) { return (i == 0); })))
        {
            continue;
        }

        const auto binaryQuals = getBinaryQuals(alignmentData.read);

        auto isBadFlank = [&](const int cstart, const int cstop) -> bool
        {
            const int span(cstop - cstart);
            if (span < minFlankLength)
            {
                return true;
            }
            const auto qiter(binaryQuals.begin() + cstart);
            const double highQBaseFraction = mean(qiter, qiter + span);
            return (highQBaseFraction < minFlankHighQBaseFraction);
        };

        // Check left flank quality
        {
            const int cstop(nodeAlignmentLengths[0]);
            const int cstart(std::max(0, cstop - minFlankLength * 2));
            if (isBadFlank(cstart, cstop))
            {
                continue;
            }
        }

        // Check right flank quality
        {
            const int cstart(nodeAlignmentLengths[0] + nodeAlignmentLengths[1]);
            const int readLength(binaryQuals.size());
            const int cstop(std::min(readLength, cstart + minFlankLength * 2));
            if (isBadFlank(cstart, cstop))
            {
                continue;
            }
        }

        numSpanningReads += 1;
    }

    return numSpanningReads;
}

/// \brief Get MotifAndPurityData from RFC1 locus reads, but attempting to exclude spanning reads
///
/// \param[in] alleleRepeatMotifCounts RFC1 repeat counts for each allele of the sample as predicted by ExpansionHunter
///
/// \param[in] expectedMotifSize: Size of the repeat motif at the target locus
///
/// \param[in] pathogenicMotifs Container of pathogenic motif strings
///
MotifAndPurityData getMotifAndPurityDataNoSpan(
    const locus::AlignmentBuffer& alignmentBuffer, const std::vector<unsigned>& alleleRepeatMotifCounts,
    const unsigned expectedMotifSize, const std::vector<std::string>& pathogenicMotifs)
{
    // If both alleles are predicted to be expanded by EH, then this is the min number of repeats the read must span to
    // be considered 'non-spanning' evidence.
    //
    // In the original proto RFC1 caller code, it was suggested this parameter could be set as a function of depth, but
    // the motivation needs to be clarified if so.
    //
    const unsigned predefinedShortRepeatMotifCount(13);

    // Get the smaller of the two repeat alleles predicted by EH:
    const unsigned minAlleleRepeatMotifCount(
        *std::min_element(alleleRepeatMotifCounts.begin(), alleleRepeatMotifCounts.end()));

    //
    const unsigned minRepeatMotifSpan = std::min(minAlleleRepeatMotifCount, predefinedShortRepeatMotifCount) + 2;

    return getMotifAndPurityData(alignmentBuffer, expectedMotifSize, minRepeatMotifSpan, pathogenicMotifs);
}

/// Loci where EH predicts at least one allele where the motif count is at least this long are treated as expanded
///
unsigned getMinExpansionRepeatMotifCount(const int readLength, const unsigned expectedMotifSize)
{
    return readLength / expectedMotifSize;
}

/// Return the total observed motif counts
///
unsigned getTotalMotifCount(const MotifObservationMap_t& motifMap)
{
    unsigned totalMotifCount(0);
    for (auto& v : motifMap | boost::adaptors::map_values)
    {
        totalMotifCount += v.count;
    }
    return totalMotifCount;
}

/// \brief Return the total fraction of pathogenic motifs
///
/// \param[in] pathogenicMotifs Container of pathogenic motif strings
///
double
getPathogenicMotifFraction(const MotifObservationMap_t& motifMap, const std::vector<std::string>& pathogenicMotifs)
{
    double totalWeightedFrac(0);
    for (const auto& motif : pathogenicMotifs)
    {
        const auto miter(motifMap.find(motif));
        if (miter != motifMap.end())
        {
            totalWeightedFrac += miter->second.weightedFrac;
        }
    }
    return totalWeightedFrac;
}

/// \brief For a sample that has already been inferred to be a double-expansion, test whether it is likely to be a
/// carrier (ie has one pathogenic and one benign expansion)
///
/// \param[in] pathogenicMotifFractionPerRead A list which for each read, contains the fraction of pathogenic motifs
/// compared to all other high-quality repeat motifs in the read.
///
/// \param[in] averageDepth Average genome depth
///
/// \return True if the sample is a carrier
///
bool isExpandedCarrier(
    const MotifObservationMap_t& motifMap, const std::vector<double>& pathogenicMotifFractionPerRead,
    const double averageDepth)
{
    if (motifMap.size() < 2)
    {
        return false;
    }

    const unsigned minReadCount(std::round(averageDepth * 0.2));
    if (pathogenicMotifFractionPerRead.size() < minReadCount)
    {
        return false;
    }

    unsigned numPathogenicReads(0);
    for (const auto pathogenicMotifFraction : pathogenicMotifFractionPerRead)
    {
        if (pathogenicMotifFraction >= 0.7)
        {
            numPathogenicReads += 1;
        }
    }

    if (numPathogenicReads > 1)
    {
        return true;
    }

    return false;
}

///
/// \param[in] numSpanningReads Count of RFC1 expansion spanning reads
///
/// \param[in] alleleRepeatMotifCounts List of predicted expansion allele lengths for RFC1
///
/// \param[in] expectedMotifSize Size of the repeat motif at the target locus
///
/// \param[in] pathogenicMotifs Container of pathogenic motif strings
///
/// \param[in] readLength
///
/// \param[in] avrerageDepth Average genome depth
///
/// \return A struct containing (1) the RFC1 call with respect ot known associations with CANVAS (2) a text
/// description elaborating more detail related to the call.
///
RFC1Status getRFC1Status(
    const std::vector<unsigned>& alleleRepeatMotifCounts, const MotifAndPurityData& mpData,
    const MotifAndPurityData& mpDataNoSpan, const unsigned numSpanningReads, const unsigned expectedMotifSize,
    const std::vector<std::string>& pathogenicMotifs, const int readLength, const double averageDepth)
{
    // The highest count of spanning reads that can still be interpreted as a sample with an expansion on both alleles
    static const unsigned maxSpanningReadsForExpansion2(1);

    assert(averageDepth >= 0.);
    const auto minNoSpanTotalMotifCount(static_cast<unsigned>(averageDepth));

    // Get initial expansion count directly from EH genotype, any repeat which extends at least to the read length is
    // counted as expanded
    unsigned expansionCount(0);
    const unsigned minExpansionCount = getMinExpansionRepeatMotifCount(readLength, expectedMotifSize);
    for (const auto alleleRepeatMotifCount : alleleRepeatMotifCounts)
    {
        if (alleleRepeatMotifCount >= minExpansionCount)
        {
            expansionCount += 1;
        }
    }

    // When there is gc bias/low-coverage, spanning reads can be used to rescue the detection of an expansion which
    // might have been missed by EH
    if (expansionCount == 0)
    {
        const unsigned noSpanTotalMotifCount(getTotalMotifCount(mpDataNoSpan.motifMap));
        if (noSpanTotalMotifCount >= minNoSpanTotalMotifCount)
        {
            // meets the revised low-coverage criteria for an expansion
            expansionCount += 1;
        }
    }
    if ((expansionCount == 1) and (numSpanningReads <= maxSpanningReadsForExpansion2))
    {
        expansionCount += 1;
    }

    // Make final classification:
    if (expansionCount == 0)
    {
        return { RFC1CallType::normal, "no expanded allele" };
    }
    else if (expansionCount == 1)
    {
        if (mpDataNoSpan.motifMap.empty())
        {
            return { RFC1CallType::normal, "expanded allele may exist but not observed" };
        }
        const double pathogenicMotifFraction(getPathogenicMotifFraction(mpDataNoSpan.motifMap, pathogenicMotifs));
        if (pathogenicMotifFraction >= 0.8)
        {
            return { RFC1CallType::carrier, "1 expanded pathogenic allele, 1 short reference allele" };
        }
        else
        {
            return { RFC1CallType::normal, "1 expanded benign allele, 1 short reference allele" };
        }
    }
    else if (expansionCount == 2)
    {
        if (mpData.motifMap.empty())
        {
            return { RFC1CallType::normal, "2 expanded alleles may exist but not observed" };
        }
        const double pathogenicMotifFraction(getPathogenicMotifFraction(mpData.motifMap, pathogenicMotifs));
        if (pathogenicMotifFraction >= 0.8)
        {
            return { RFC1CallType::affected, "2 expanded pathogenic alleles" };
        }
        else
        {
            // first check if 1 expanded pathogenic motif allele and 1 expanded benign motif allele
            if (isExpandedCarrier(mpData.motifMap, mpData.pathogenicMotifFractionPerRead, averageDepth))
            {
                return { RFC1CallType::carrier, "1 expanded pathogenic allele, 1 expanded benign allele" };
            }

            // if not, do some general classification
            if (pathogenicMotifFraction >= 0.3)
            {
                return { RFC1CallType::potential_carrier, "2 expanded alleles with >30% pathogenic kmers" };
            }
            else
            {
                return { RFC1CallType::normal, "2 expanded alleles (possibly reference)" };
            }
        }
    }
    else
    {
        const std::string message("Illegal expansion count: " + std::to_string(expansionCount));
        throw std::logic_error(message);
    }
}

}

void runRFC1MotifAnalysis(const locus::AlignmentBuffer& alignmentBuffer, LocusFindings& locusFindings)
{
    // Note that the 'use_rotation' and 'use_spanning' options from the proto version of this method are both fixed to
    // true here.
    //

    // Hard coded parameters for RFC1 locus:
    static const unsigned expectedMotifSize(5);
    static const std::vector<std::string> pathogenicMotifs { "AAGGG", "ACAGG" };

    // RFC1 motif analysis loci are constrained to only one variant (this is enforced with an error message when loading
    // the catalog)
    assert(locusFindings.findingsForEachVariant.size() == 1);

    const double averageDepth(locusFindings.stats.depth());
    const int readLength(locusFindings.stats.meanReadLength());

    // Extract standard EH results from repeat findings, and add RFC1 results back in here as a final step:
    RepeatFindings& repeatFindings(
        *dynamic_cast<RepeatFindings*>(locusFindings.findingsForEachVariant.begin()->second.get()));

    // Get standard motif map (including any spanning reads)
    //
    static const unsigned standardMinRepeatMotifSpan(0);
    const auto mpData
        = getMotifAndPurityData(alignmentBuffer, expectedMotifSize, standardMinRepeatMotifSpan, pathogenicMotifs);

    // Get 'no-spanning' motif map (attempting to exclude spanning reads)
    //
    const std::vector<unsigned> alleleRepeatMotifCounts(getAlleleRepeatMotifCounts(repeatFindings));
    const auto mpDataNoSpan
        = getMotifAndPurityDataNoSpan(alignmentBuffer, alleleRepeatMotifCounts, expectedMotifSize, pathogenicMotifs);

    const unsigned numSpanningReads(countSpanningReads(alignmentBuffer));
    const RFC1Status rfc1Status = getRFC1Status(
        alleleRepeatMotifCounts, mpData, mpDataNoSpan, numSpanningReads, expectedMotifSize, pathogenicMotifs,
        readLength, averageDepth);

    // Report analysis results out to the appropriate RepeatFindings structure:
    repeatFindings.setRFC1Status(rfc1Status);
}

}
