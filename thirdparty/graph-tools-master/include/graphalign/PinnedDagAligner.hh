//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Roman Petrovski <RPetrovski@illumina.com>
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

#pragma once

#include <list>
#include <map>
#include <string>

#include <boost/throw_exception.hpp>

#include "graphalign/DagAlignerAffine.hh"
#include "graphalign/GraphAlignment.hh"
#include "graphalign/Operation.hh"
#include "graphalign/dagAligner/BaseMatchingPenaltyMatrix.hh"
#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

namespace graphtools
{

using PathAndAlignment = std::pair<Path, Alignment>;

template <bool penalizeMove, bool clipFront = true>
class BaseMatchingDagAligner : public graphalign::dagAligner::Aligner<
                                   graphalign::dagAligner::AffineAlignMatrixVectorized<
                                       graphalign::dagAligner::BaseMatchingPenaltyMatrix, penalizeMove>,
                                   clipFront>
{
    typedef graphalign::dagAligner::BaseMatchingPenaltyMatrix PenaltyMatrix;
    typedef graphalign::dagAligner::Score Score;

public:
    BaseMatchingDagAligner(const PenaltyMatrix& penaltyMatrix, Score gapOpen, Score gapExt)
        : graphalign::dagAligner::Aligner<
              graphalign::dagAligner::AffineAlignMatrixVectorized<PenaltyMatrix, penalizeMove>, clipFront>(
              penaltyMatrix, gapOpen, gapExt)
    {
    }

    BaseMatchingDagAligner(graphalign::dagAligner::Score match, Score mismatch, Score gapOpen, Score gapExt)
        : graphalign::dagAligner::Aligner<
              graphalign::dagAligner::AffineAlignMatrixVectorized<PenaltyMatrix, penalizeMove>, clipFront>(
              PenaltyMatrix(match, mismatch), gapOpen, gapExt)
    {
    }
};

/**
 * Performs alignment of query pieces that start or end at the seed in the graph
 */
class PinnedDagAligner
{
    typedef std::pair<int, int> Edge;
    typedef std::vector<Edge> Edges;
    typedef graphalign::dagAligner::Cigar Cigar;
    BaseMatchingDagAligner<true, false> aligner_;

    static void appendOperation(OperationType type, uint32_t length, std::list<Operation>& operations)
    {
        if (operations.empty() || operations.back().type() != type)
        {
            operations.push_back(Operation(type, length));
        }
        else
        {
            operations.back() = Operation(type, operations.back().length() + length);
        }
    }

    template <typename GraphT, typename PathT>
    static void parseGraphCigar(
        const GraphT& graph, const graphalign::dagAligner::Cigar& cigar, PathT& path, std::list<Operation>& operations)
    {
        using namespace graphalign::dagAligner;
        for (const Cigar::Operation& op : cigar)
        {
            if (Cigar::NODE_START == op.code_)
            {
                if (path.lastNodeId() != op.value_
                    || std::size_t(path.endPosition()) == graph.nodeSeq(op.value_).length())
                {
                    path.extendEndToNode(op.value_);
                }
            }
            else if (Cigar::NODE_END == op.code_)
            {
            }
            else
                switch (op.code_)
                {
                case Cigar::MATCH:
                {
                    appendOperation(OperationType::kMatch, op.value_, operations);
                    path.shiftEndAlongNode(op.value_);
                    break;
                }
                case Cigar::MISMATCH:
                {
                    appendOperation(OperationType::kMismatch, op.value_, operations);
                    path.shiftEndAlongNode(op.value_);
                    break;
                }
                case Cigar::INSERT:
                {
                    appendOperation(OperationType::kInsertionToRef, op.value_, operations);
                    break;
                }
                case Cigar::SOFT_CLIP:
                {
                    appendOperation(OperationType::kSoftclip, op.value_, operations);
                    break;
                }
                case Cigar::DELETE:
                {
                    appendOperation(OperationType::kDeletionFromRef, op.value_, operations);
                    path.shiftEndAlongNode(op.value_);
                    break;
                }
                default:
                {
                    throw std::logic_error(std::string("Unexpected graph cigar operation:") + std::to_string(op.code_));
                }
                }
        }
    }

    typedef int MappedId;

    void unmapNodeIds(const std::map<MappedId, NodeId>& originalIds, graphalign::dagAligner::Cigar& cigar)
    {
        using namespace graphalign::dagAligner;
        for (Cigar::Operation& op : cigar)
        {
            if (Cigar::NODE_START == op.code_ || Cigar::NODE_END == op.code_)
            {
                op.value_ = originalIds.at(op.value_);
            }
        }
    }

    /**
     * \brief Tests if the cigar first node is a repeat expansion and corrects
     *        cigar to ensure that fraction of the first expansion is interpreted correctly
     */
    template <typename PathT>
    void fixFirstNodeExpansion(
        const std::vector<MappedId>& nodeIds, const std::map<MappedId, NodeId>& originalIds, const PathT& seedPath,
        graphalign::dagAligner::Cigar& cigar)
    {
        using namespace graphalign::dagAligner;

        if (seedPath.lastNodeId() != originalIds.at(cigar.firstNode()) || MappedId(cigar.firstNode()) != nodeIds[0])
        {
            const int skipLen
                = seedPath.graphRawPtr()->nodeSeq(seedPath.lastNodeId()).length() - seedPath.endPosition();
            if (0 > skipLen)
            {
                throw std::logic_error("invalid distance to skip");
            }
            if (skipLen)
            {
                const Cigar::Operation skipFirstNode[]
                    = { Cigar::Operation(
                            Cigar::NODE_START, nodeIds[0]), // endPosition is on the base that belongs to the path...
                        Cigar::Operation(Cigar::DELETE, skipLen), Cigar::Operation(Cigar::NODE_END, nodeIds[0]) };

                cigar.insert(
                    cigar.begin(), skipFirstNode, skipFirstNode + sizeof(skipFirstNode) / sizeof(skipFirstNode[0]));
            }
        }
    }

public:
    explicit PinnedDagAligner(
        const int32_t matchScore, const int32_t mismatchScore, const int32_t gapOpenScore, const int32_t gapExtendScore)
        : aligner_(matchScore, mismatchScore, gapOpenScore, gapExtendScore)
    {
    }

    std::list<PathAndAlignment>
    prefixAlign(const Path& seedPath, const std::string& queryPiece, size_t extensionLen, int& score)
    {
        using namespace graphalign::dagAligner;
        std::vector<MappedId> nodeIds;
        Edges edges;
        std::string target;

        // when repeat expansions are unrolled each copy gets a unique id, so, all
        // ids have to be remapped
        std::map<MappedId, NodeId> originalIds;

        bfsDiscoverEdges(
            *seedPath.graphRawPtr(), seedPath.nodeIds().back(), seedPath.endPosition(), extensionLen, nodeIds, edges,
            target, originalIds);
        edges.push_back(Edge(target.length(), target.length()));

        std::list<PathAndAlignment> ret;
        if (!target.empty())
        {
            const EdgeMap alignerEdges(edges, nodeIds);

            aligner_.align(queryPiece.begin(), queryPiece.end(), target.begin(), target.end(), alignerEdges);

            std::vector<Cigar> cigars;
            Score secondBestScore = 0;
            Score bestScore = aligner_.backtrackAllPaths<false>(alignerEdges, cigars, secondBestScore);

            score = bestScore;
            for (Cigar& cigar : cigars)
            {
                fixFirstNodeExpansion(nodeIds, originalIds, seedPath, cigar);

                unmapNodeIds(originalIds, cigar);

                Path path = seedPath;
                std::list<Operation> operations;
                parseGraphCigar(*seedPath.graphRawPtr(), cigar, path, operations);

                ret.push_back(PathAndAlignment(path, Alignment(seedPath.seq().length(), operations)));
            }
        }

        return ret;
    }

    std::list<PathAndAlignment>
    suffixAlign(const Path& seedPath, std::string queryPiece, size_t extensionLen, int& score)
    {
        using namespace graphalign::dagAligner;
        std::vector<MappedId> nodeIds;
        Edges edges;
        std::string target;

        // when repeat expansions are unrolled each copy gets a unique id, so, all
        // ids have to be remapped
        std::map<MappedId, NodeId> originalIds;

        ReverseGraph rg(*seedPath.graphRawPtr());
        bfsDiscoverEdges(
            rg, seedPath.nodeIds().front(),
            // endPosition is on the base that belongs to the path...
            ConstReversePath(seedPath).endPosition(), extensionLen, nodeIds, edges, target, originalIds);
        edges.push_back(Edge(target.length(), target.length()));

        std::list<PathAndAlignment> ret;
        if (!target.empty())
        {
            const EdgeMap alignerEdges(edges, nodeIds);

            std::reverse(queryPiece.begin(), queryPiece.end());
            aligner_.align(queryPiece.begin(), queryPiece.end(), target.begin(), target.end(), alignerEdges);

            std::vector<Cigar> cigars;
            Score secondBestScore = 0;
            Score bestScore = aligner_.backtrackAllPaths<false>(alignerEdges, cigars, secondBestScore);

            score = bestScore;
            for (Cigar& cigar : cigars)
            {
                fixFirstNodeExpansion(nodeIds, originalIds, ConstReversePath(seedPath), cigar);

                unmapNodeIds(originalIds, cigar);

                Path path = seedPath;
                ReversePath rp(path);
                std::list<Operation> operations;
                parseGraphCigar(rg, cigar, rp, operations);
                operations.reverse();

                // reversed alignments always start at the beginning of the path because
                // the seed path gets start-extended to incorporate them
                ret.push_back(PathAndAlignment(path, Alignment(0, operations)));
            }
        }

        return ret;
    }

private:
    template <typename GraphT>
    static std::map<NodeId, int> extractSubgraph(
        const GraphT& graph, const NodeId startNodeId, const std::size_t startNodeOffset, const std::size_t seqLen)
    {
        std::map<NodeId, int> nodeStartSeqOffset;

        if (graph.nodeSeq(startNodeId).length() == startNodeOffset)
        {
            // flag empty start node in a special way
            nodeStartSeqOffset[startNodeId] = -1;
        }
        else
        {
            // first node start is at the start of the sequence
            nodeStartSeqOffset[startNodeId] = 0;
        }

        // nodes to be visited by bfs
        std::deque<NodeId> shouldVisit(1, startNodeId);

        // extract longest subgraph of nodes such that each node begins within
        // seqLen from the start of the start node
        while (!shouldVisit.empty())
        {
            const NodeId currentNodeId = shouldVisit.front();
            shouldVisit.pop_front();

            const std::string& currentNodeSeq = graph.nodeSeq(currentNodeId);
            const int currentNodeSeqOffset = nodeStartSeqOffset[currentNodeId];

            // avoid dealing with individual node startNodeOffset (only start node has it)
            // by pretending the sequence starts at the start node start
            if (seqLen + startNodeOffset - std::max(0, currentNodeSeqOffset) > currentNodeSeq.length())
            {
                const int successorSeqOffset
                    = -1 == currentNodeSeqOffset ? 0 : currentNodeSeqOffset + currentNodeSeq.length();
                // sequence does not terminate at this node, enqueue successors
                const std::set<NodeId>& successors = graph.successors(currentNodeId);
                for (const NodeId successorId : successors)
                {
                    const auto seenSuccessor = nodeStartSeqOffset.find(successorId);
                    if (nodeStartSeqOffset.end() == seenSuccessor || seenSuccessor->second > successorSeqOffset)
                    {
                        // some successors will end up listed in shouldVisit more than once at a time
                        shouldVisit.push_back(successorId);
                        nodeStartSeqOffset[successorId] = successorSeqOffset;
                    }
                }
            }
        }

        return nodeStartSeqOffset;
    }

    typedef std::pair<MappedId, MappedId> IdEdge;
    typedef std::vector<IdEdge> IdEdges;

    /**
     * \brief expands repeats up to remainder of sequnce length
     * \return pairs of mapped node ids indicating an edge between them.
     * \postcondition result array is ordered by successor id then by predecessor id
     */
    template <typename GraphT>
    static IdEdges unrollRepeats(
        const GraphT& graph, const std::size_t seqLen, const std::map<NodeId, int>& nodeStartSeqOffset,
        std::map<MappedId, NodeId>& originalIds, std::multimap<NodeId, MappedId>& mappedIds)
    {
        IdEdges idEdges;
        for (const auto& nodeIdOffset : nodeStartSeqOffset)
        {
            if (int(seqLen) <= nodeIdOffset.second)
            {
                throw std::logic_error("node should not be in the subgraph");
            }

            const std::set<NodeId>& successors = graph.successors(nodeIdOffset.first);
            if (successors.end() != successors.find(nodeIdOffset.first))
            {
                const std::size_t nodeSeqLen = graph.nodeSeq(nodeIdOffset.first).length();
                // don't forget empty repeat start node has special sequence offset -1
                for (std::size_t lenLeft = seqLen - std::max(0, nodeIdOffset.second); lenLeft;)
                {
                    // chain unrolled repeat nodes together
                    IdEdge edge(mappedIds.size(), mappedIds.size() + 1);
                    idEdges.push_back(edge);
                    mappedIds.emplace(nodeIdOffset.first, originalIds.size());
                    originalIds.emplace(originalIds.size(), nodeIdOffset.first);
                    lenLeft -= std::min(lenLeft, nodeSeqLen);
                }

                // since edges point forward, above loop always produces one more edge than we need
                idEdges.pop_back();
            }
            else
            {
                mappedIds.emplace(nodeIdOffset.first, originalIds.size());
                originalIds.emplace(originalIds.size(), nodeIdOffset.first);
            }
        }

        linkPredecessors(graph, originalIds, mappedIds, idEdges);

        // group edges by successor node
        std::sort(idEdges.begin(), idEdges.end(), [](const IdEdge& left, const IdEdge& right) {
            return left.second < right.second || (left.second == right.second && left.first < right.first);
        });

        return idEdges;
    }

    static void dfsExtractOrderedNodeIds(
        const MappedId currentId, const IdEdges& idEdges, const std::vector<std::size_t>& idEdgesIndex,
        std::vector<char>& seenNodes, std::vector<MappedId>& nodeIds)
    {
        if (!seenNodes.at(currentId))
        {
            seenNodes[currentId] = true;
            for (std::size_t predOffset = idEdgesIndex.at(currentId); idEdgesIndex.at(currentId + 1) != predOffset;
                 ++predOffset)
            {
                if (idEdges.at(predOffset).second != currentId)
                {
                    throw std::logic_error(
                        "dfsExtractOrderedNodeIds: Invalid edge for node " + std::to_string(currentId));
                }

                dfsExtractOrderedNodeIds(idEdges[predOffset].first, idEdges, idEdgesIndex, seenNodes, nodeIds);
            }
            nodeIds.push_back(currentId);
        }
    }
    /**
     * \brief dfs in order to produce the topological ordering with start node on top
     */
    static std::vector<MappedId>
    extractOrderedNodeIds(const IdEdges& idEdges, const std::vector<std::size_t>& idEdgesIndex)
    {
        std::vector<MappedId> nodeIds;
        std::vector<char> seenNodes(idEdgesIndex.size(), false);
        for (std::size_t mappedId = 0; idEdgesIndex.size() - 1 != mappedId; ++mappedId)
        {
            dfsExtractOrderedNodeIds(mappedId, idEdges, idEdgesIndex, seenNodes, nodeIds);
        }

        return nodeIds;
    }

    /**
     * \brief self-repeat edges already in idEdges, add all mapped predecessor edges to each first expansion
     *        and non-repeat nodes
     */
    template <typename GraphT>
    static void linkPredecessors(
        const GraphT& graph, const std::map<MappedId, NodeId>& originalIds,
        const std::multimap<NodeId, MappedId>& mappedIds, IdEdges& idEdges)
    {
        for (MappedId mappedId = 0; MappedId(mappedIds.size()) != mappedId;)
        {
            const NodeId originalId = originalIds.at(mappedId);

            bool selfRepeat = false;
            for (const NodeId& predecessorId : graph.predecessors(originalId))
            {
                if (predecessorId != originalId)
                {
                    auto mappedPredIds = mappedIds.equal_range(predecessorId);

                    // other nodes, insert edges for each predecessor instance
                    // empty ranges indicate predecessors that are not part of subgraph
                    while (mappedPredIds.second != mappedPredIds.first)
                    {
                        IdEdge edge(mappedPredIds.first->second, mappedId);
                        idEdges.push_back(edge);
                        ++mappedPredIds.first;
                    }
                }
                else
                {
                    selfRepeat = true;
                }
            }

            if (selfRepeat)
            {
                // skip all instances of self repeat so that only first node gets edges
                // from predecessors
                auto mappedIdsRange = mappedIds.equal_range(originalId);
                mappedId += std::distance(mappedIdsRange.first, mappedIdsRange.second);
            }
            else
            {
                ++mappedId;
            }
        }
    }

    static std::vector<std::size_t> indexEdges(const IdEdges& idEdges, const std::map<MappedId, NodeId>& originalIds)
    {
        std::vector<std::size_t> index;
        index.reserve(idEdges.size() + 1);
        index.push_back(0);
        if (!idEdges.size())
        {
            // no edges. Single-node graph, just close the index
            index.push_back(0);
        }
        else
        {
            MappedId lastId = -1;
            for (const IdEdge& edge : idEdges)
            {
                while (lastId != edge.second)
                {
                    index.push_back(index.back());
                    ++lastId;
                }
                ++index.back();
            }
            // add nodes without predecessors to index
            while (lastId != MappedId(originalIds.rbegin()->first))
            {
                index.push_back(index.back());
                ++lastId;
            }
        }

        return index;
    }

    template <typename GraphT>
    static std::string buildTargetSequence(
        const GraphT& graph, const std::size_t startNodeOffset, const std::vector<MappedId>& nodeIds,
        const std::map<MappedId, NodeId>& originalIds, const std::map<NodeId, int>& nodeStartSeqOffset,
        const std::vector<std::size_t>& idEdgesIndex, const IdEdges& idEdges, std::vector<Edge>& edges)
    {
        std::string target;
        std::vector<int> mappedIdEndOffset(nodeIds.size(), 0);
        // when first node is a repeat expansion fully consumed by seed, just pretend that query start
        // at the beginning of the node rather than after the end...
        std::size_t startOffset = startNodeOffset;
        for (const MappedId mappedId : nodeIds)
        {
            const NodeId originalId = originalIds.at(mappedId);
            const std::string& nodeSeq = graph.nodeSeq(originalId);
            if (!startOffset || -1 != nodeStartSeqOffset.at(originalId))
            {
                const std::size_t nodePartLen = nodeSeq.length() - startOffset;
                if (!nodePartLen)
                {
                    throw std::logic_error("Empty node in expanded subgraph and it is not the first one");
                }

                for (std::size_t predOffset = idEdgesIndex[mappedId]; idEdgesIndex[mappedId + 1] != predOffset;
                     ++predOffset)
                {
                    if (idEdges[predOffset].second != mappedId)
                    {
                        throw std::logic_error("bfsDiscoverEdges: Invalid edge");
                    }
                    Edge edge(mappedIdEndOffset.at(idEdges[predOffset].first), target.length());
                    edges.push_back(edge);
                }
                target += nodeSeq.substr(startOffset, nodePartLen);
            }
            mappedIdEndOffset[mappedId] = target.length() - 1;
            startOffset = 0;
        }

        return target;
    }

    /**
     * \param seqLen    extensionLen + offset in the first node. This simplifies
     * the process
     *                  by assuming that sequence starts at the node start
     */
    template <typename GraphT>
    static void bfsDiscoverEdges(
        const GraphT& graph, const NodeId startNodeId, const std::size_t startNodeOffset, const std::size_t seqLen,
        std::vector<MappedId>& nodeIds, std::vector<Edge>& edges, std::string& target,
        std::map<MappedId, NodeId>& originalIds)
    {
        // length of shortest path to the first node character
        const std::map<NodeId, int> nodeStartSeqOffset = extractSubgraph(graph, startNodeId, startNodeOffset, seqLen);

        // Repeat expansions need to be unrolled, so we create unique id for each unrolled instance and map
        // them to original ids
        std::multimap<NodeId, MappedId> mappedIds;
        const IdEdges idEdges
            = unrollRepeats(graph, startNodeOffset + seqLen, nodeStartSeqOffset, originalIds, mappedIds);

        // index the edge array so that for each mappedId there is a pair of entries index[mappedId],
        // index[mappedId+1] which contains offsets of all predecessor edges in the idEdges;
        const std::vector<std::size_t> idEdgesIndex = indexEdges(idEdges, originalIds);

        nodeIds = extractOrderedNodeIds(idEdges, idEdgesIndex);
        if (nodeIds.size() != mappedIds.size())
        {
            throw std::logic_error(
                "Invalid number of nodeIds entries " + std::to_string(nodeIds.size()) + " expected "
                + std::to_string(mappedIds.size()));
        }

        if (originalIds.at(nodeIds.at(0)) != startNodeId)
        {
            // nodeIds must be topologically sorted. Since we're discovering our subgraph from the
            // start node, its first expansion must sort to the top.
            throw std::logic_error("First expansion of start node must be the first node");
        }

        // extract target sequence in the proper order
        target = buildTargetSequence(
            graph, startNodeOffset, nodeIds, originalIds, nodeStartSeqOffset, idEdgesIndex, idEdges, edges);

        if (-1 == nodeStartSeqOffset.at(startNodeId))
        {
            // if the first node is an empty repeat expansion (consumed by seed), remove it as it has no
            // corresponding sequence in the target.
            nodeIds.erase(nodeIds.begin());
        }
    }
};
}
