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

#include "region_analysis/RegionAnalyzer.hh"

#include <cassert>
#include <iostream>

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphutils/SequenceOperations.hh"

#include "alignment/AlignmentFilters.hh"
#include "alignment/AlignmentTweakers.hh"
#include "alignment/GraphAlignmentOperations.hh"

using graphtools::Operation;
using graphtools::OperationType;
using graphtools::splitStringByDelimiter;
using reads::LinearAlignmentStats;
using reads::Read;
using std::list;
using std::string;
using std::vector;

void outputAlignments(const list<GraphAlignment>& alignments)
{
    for (const auto& alignment : alignments)
    {
        std::cerr << "\t\t" << alignment << std::endl;
    }
}

static string indentMultilineString(const string str, int32_t indentation_len)
{
    string indented_str;
    const vector<string> lines = splitStringByDelimiter(str, '\n');
    for (auto& line : lines)
    {
        if (!indented_str.empty())
        {
            indented_str += '\n';
        }
        indented_str += string(indentation_len, ' ') + line;
    }

    return indented_str;
}

static void outputAlignedRead(const Read& read, boost::optional<GraphAlignment> alignment, std::ostream& out)
{
    const int32_t indentation_size = 2;
    const string spacer(indentation_size, ' ');
    out << spacer << "- name: " << read.read_id << std::endl;
    if (alignment)
    {
        out << spacer << "  path: " << alignment->path() << std::endl;
        out << spacer << "  graph_cigar: " << alignment->generateCigar() << std::endl;
        out << spacer << "  alignment: |" << std::endl;
        const string alignment_encoding = prettyPrint(*alignment, read.sequence);
        out << indentMultilineString(alignment_encoding, 3 * indentation_size) << std::endl;
    }
}

void RegionAnalyzer::processMates(reads::Read read1, reads::Read read2)
{
    boost::optional<GraphAlignment> read1Alignment = alignRead(read1);
    outputAlignedRead(read1, read1Alignment, alignmentStream_);

    boost::optional<GraphAlignment> read2Alignment = alignRead(read2);
    outputAlignedRead(read2, read2Alignment, alignmentStream_);

    for (auto& repeatAnalyzer : repeatAnalyzers_)
    {
        repeatAnalyzer.processMates(
            read1.readId(), read1.sequence, read1Alignment, read2.readId(), read2.sequence, read2Alignment);
    }
}

boost::optional<GraphAlignment> RegionAnalyzer::alignRead(Read& read) const
{
    OrientationPrediction predictedOrientation = orientationPredictor_.predict(read.sequence);

    if (predictedOrientation == OrientationPrediction::kAlignsInReverseComplementOrientation)
    {
        read.sequence = graphtools::reverseComplement(read.sequence);
    }
    else if (predictedOrientation == OrientationPrediction::kDoesNotAlign)
    {
        return boost::optional<GraphAlignment>();
    }

    const list<GraphAlignment> alignments = graphAligner_.align(read.sequence);

    if (alignments.empty())
    {
        return boost::optional<GraphAlignment>();
    }

    GraphAlignment canonicalAlignment = computeCanonicalAlignment(alignments);

    if (checkIfPassesAlignmentFilters(canonicalAlignment))
    {
        const int kShrinkLength = 10;
        shrinkUncertainPrefix(kShrinkLength, read.sequence, canonicalAlignment);
        shrinkUncertainSuffix(kShrinkLength, read.sequence, canonicalAlignment);

        return canonicalAlignment;
    }
    else
    {
        return boost::optional<GraphAlignment>();
    }
}

bool RegionAnalyzer::checkIfPassesAlignmentFilters(const GraphAlignment& alignment) const
{
    const Operation& firstOperation = alignment.alignments().front().operations().front();
    const int frontSoftclipLen = firstOperation.type() == OperationType::kSoftclip ? firstOperation.queryLength() : 0;

    const Operation& lastOperation = alignment.alignments().back().operations().back();
    const int backSoftclipLen = lastOperation.type() == OperationType::kSoftclip ? lastOperation.queryLength() : 0;

    const int clippedQueryLength = alignment.queryLength() - frontSoftclipLen - backSoftclipLen;
    const int referenceLength = alignment.referenceLength();

    const int percentQueryMatches = (100 * alignment.numMatches()) / clippedQueryLength;
    const int percentReferenceMatches = (100 * alignment.numMatches()) / referenceLength;

    if (percentQueryMatches >= 80 && percentReferenceMatches >= 80)
    {
        return true;
    }
    else
    {
        return false;
    }
}

RegionFindings RegionAnalyzer::genotype()
{
    RegionFindings regionResults;

    for (auto& repeatAnalyzer : repeatAnalyzers_)
    {
        RepeatFindings repeatFindings = repeatAnalyzer.analyzeCommonRepeat();
        regionResults.emplace(std::make_pair(repeatAnalyzer.repeatId(), repeatFindings));
    }

    return regionResults;
}

vector<std::unique_ptr<RegionAnalyzer>> initializeRegionAnalyzers(
    const RegionCatalog& RegionCatalog, double haplotypeDepth, int readLength, const string& alignerName,
    std::ostream& alignmentStream)
{
    vector<std::unique_ptr<RegionAnalyzer>> regionAnalyzers;

    for (const auto& regionIdAndRegionSpec : RegionCatalog)
    {

        const RegionSpec& regionSpec = regionIdAndRegionSpec.second;

        GraphAlignmentHeuristicsParameters alignmentParams;

        // alignmentStream << regionSpec.regionId() << ":" << std::endl;
        regionAnalyzers.emplace_back(
            new RegionAnalyzer(regionSpec, haplotypeDepth, readLength, alignmentStream, alignerName, alignmentParams));
    }

    return regionAnalyzers;
}