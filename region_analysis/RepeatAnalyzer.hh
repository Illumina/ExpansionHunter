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

#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

#include "classification/alignment_classifier.h"
#include "common/common.h"
#include "genotyping/RepeatGenotype.hh"
#include "reads/read.h"
#include "region_analysis/RepeatFindings.hh"

class RepeatAnalyzer
{
public:
    RepeatAnalyzer(
        const std::string& repeatId, graphtools::Graph graph, AlleleCount expectedAlleleCount,
        graphtools::NodeId repeatNodeId, double haplotypeDepth, int32_t maxNumUnitsInRead)
        : repeatId_(repeatId)
        , expectedAlleleCount_(expectedAlleleCount)
        , repeatNodeId_(repeatNodeId)
        , maxNumUnitsInRead_(maxNumUnitsInRead)
        , haplotypeDepth_(haplotypeDepth)
        , repeatUnit_(graph.nodeSeq(repeatNodeId_))
        , alignmentClassifier_(graph, repeatNodeId)
    {
    }

    const std::string& repeatId() const { return repeatId_; }
    void processMates(
        const std::string& readName, const std::string& readSequence,
        boost::optional<graphtools::GraphAlignment> readAlignment, const std::string& mateName,
        const std::string& mateSequence, boost::optional<graphtools::GraphAlignment> mateAlignment);

    RepeatFindings analyzeCommonRepeat() const;

    bool operator==(const RepeatAnalyzer& other) const;

private:
    boost::optional<RepeatGenotype> genotypeCommonRepeat(const std::vector<int32_t>& candidateAlleleSizes) const;
    reads::RepeatAlignmentStats classifyReadAlignment(const graphtools::GraphAlignment& alignment);
    void summarizeAlignmentsToReadCounts(const reads::RepeatAlignmentStats& repeatAlignmentStats);
    std::vector<int32_t> generateCandidateAlleleSizes() const;

    const std::string repeatId_;
    const AlleleCount expectedAlleleCount_;
    const graphtools::NodeId repeatNodeId_;
    const int32_t maxNumUnitsInRead_;
    const double haplotypeDepth_;

    const std::string repeatUnit_;
    RepeatAlignmentClassifier alignmentClassifier_;

    CountTable countsOfSpanningReads_;
    CountTable countsOfFlankingReads_;
    CountTable countsOfRepeatReads_;
};
