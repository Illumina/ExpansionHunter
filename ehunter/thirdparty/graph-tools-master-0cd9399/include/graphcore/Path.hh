//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
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

#pragma once

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

#include "graphcore/Graph.hh"

namespace graphtools
{

/**
 * A path in a sequence graph is given by (1) a sequence of nodes and (2) start/end position on the first/last node. The
 * start/end positions are 0-based and form a half-open interval.
 */
class Path
{
public:
    typedef std::vector<NodeId>::const_iterator const_iterator;

    /// Checks if the inputs define a well-formed path.
    Path(const Graph* graph_raw_ptr, int32_t start_position, const std::vector<NodeId>& nodes, int32_t end_position);
    ~Path() = default;

    bool operator==(const Path& other) const;
    bool operator<(const Path& path) const;

    const_iterator begin() const { return nodes_.begin(); }
    const_iterator end() const { return nodes_.end(); }

    /// Ids of nodes overlapped by the path
    const std::vector<NodeId>& nodeIds() const { return nodes_; }
    size_t numNodes() const { return nodeIds().size(); }

    /// Sequence of the entire path
    std::string seq() const;

    /// Piece of node sequence that the path overlaps
    std::string getNodeSeq(size_t node_index) const;
    const Graph* graphRawPtr() const { return graph_raw_ptr_; }
    std::string encode() const;
    int32_t startPosition() const { return start_position_; }
    int32_t endPosition() const { return end_position_; }
    size_t length() const;
    NodeId getNodeIdByIndex(const size_t node_index) const { return nodes_[node_index]; }
    int32_t getStartPositionOnNodeByIndex(size_t node_index) const;
    int32_t getEndPositionOnNodeByIndex(size_t node_index) const;
    size_t getNodeOverlapLengthByIndex(size_t node_index) const;

    bool checkOverlapWithNode(NodeId node_id) const;
    int32_t getDistanceFromPathStart(NodeId node, int32_t offset) const;

    // Path modifiers
    // Moves start position by a specified length; the path gets longer/shorter if shift_len is positive/negative
    // respectively.
    void shiftStartAlongNode(int32_t shift_len);
    // Moves end position by a specified length; the path gets longer/shorter if shift_len is positive/negative
    // respectively.
    void shiftEndAlongNode(int32_t shift_len);
    // Moves path start to the end of the specified node. The new node must be a predecessor of the previous start node.
    void extendStartToNode(NodeId node_id);
    // Moves path start to the start of the specified node. The new node must be a predecessor of the previous start
    // node.
    void extendStartToIncludeNode(NodeId node_id);
    // Moves path end to the start of the specified node. The new node must be a successor of the previous start node.
    void extendEndToNode(NodeId node_id);
    // Moves path end to the end of the specified node. The new node must be a successor of the previous start node.
    void extendEndToIncludeNode(NodeId node_id);
    // Moves path start to the start of the next node in the path.
    void removeStartNode();
    // Moves path end to the end of the previous node in the path.
    void removeEndNode();
    // Moves path start to the start of the next node if the start has zero-length overlap with the corresponding node;
    // does nothing if path spans only one node.
    void removeZeroLengthStart();
    // Moves path end to the end of the end of the previous node if the end of the path has zero-length overlap with the
    // corresponding node; does nothing if path spans only one node.
    void removeZeroLengthEnd();
    // Shortens the start of the path by a specified length.
    void shrinkStartBy(int32_t shrink_len);
    // Shortens the end of the path by a specified length.
    void shrinkEndBy(int32_t shrink_len);
    // Shortens the path by the specified lengths from each direction.
    void shrinkBy(int32_t start_shrink_len, int32_t end_shrink_len);

    NodeId firstNodeId() const { return nodeIds().front(); }
    NodeId lastNodeId() const { return nodeIds().back(); }

private:
    void assertValidity() const;
    bool isNodePositionValid(NodeId node_id, int32_t position) const;
    void assertPositionsOrdered() const;
    void assertNonEmpty() const;
    void assertFirstNodePosValid() const;
    void assertLastNodePosValid() const;
    void assertConnected() const;
    void assertThatIndexIsValid(int32_t node_index) const;

    const Graph* graph_raw_ptr_;
    int32_t start_position_;
    int32_t end_position_;
    std::vector<NodeId> nodes_;
};

std::ostream& operator<<(std::ostream& os, const Path& path);

class ReversePath
{
    Path& path_;

public:
    explicit ReversePath(Path& path)
        : path_(path)
    {
    }

    // TODO: add methods as needed

    NodeId firstNodeId() const { return path_.lastNodeId(); }
    NodeId lastNodeId() const { return path_.firstNodeId(); }

    int32_t startPosition() const
    {
        return path_.graphRawPtr()->nodeSeq(path_.lastNodeId()).length() - path_.endPosition();
    }

    int32_t endPosition() const
    {
        return path_.graphRawPtr()->nodeSeq(path_.firstNodeId()).length() - path_.startPosition();
    }

    std::string seq() const
    {
        std::string ret = path_.seq();
        std::reverse(ret.begin(), ret.end());
        return ret;
    }

    void shiftEndAlongNode(int32_t shift_len) { path_.shiftStartAlongNode(shift_len); }
    void extendEndToNode(NodeId node_id) { path_.extendStartToNode(node_id); }

    friend std::ostream& operator<<(std::ostream& os, const ReversePath& path);
};

class ConstReversePath
{
    const Path& path_;

public:
    explicit ConstReversePath(const Path& path)
        : path_(path)
    {
    }

    // TODO: add methods as needed

    NodeId firstNodeId() const { return path_.lastNodeId(); }
    NodeId lastNodeId() const { return path_.firstNodeId(); }

    int32_t startPosition() const
    {
        return path_.graphRawPtr()->nodeSeq(path_.lastNodeId()).length() - path_.endPosition();
    }

    int32_t endPosition() const
    {
        return path_.graphRawPtr()->nodeSeq(path_.firstNodeId()).length() - path_.startPosition();
    }

    const Graph* graphRawPtr() const { return path_.graphRawPtr(); }
};
}
