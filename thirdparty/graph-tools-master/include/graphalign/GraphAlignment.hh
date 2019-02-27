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

#include <initializer_list>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <regex>
#include <string>
#include <utility>
#include <vector>

#include "graphalign/LinearAlignment.hh"
#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

namespace graphtools
{
/**
 * Represents an alignment of a sequence to a graph. Graph alignment consists of a path and linear alignments for each
 * node of the path.
 */
class GraphAlignment
{
public:
    typedef size_t size_type;
    typedef std::vector<Alignment> NodeAlignments;
    typedef NodeAlignments::const_iterator const_iterator;
    GraphAlignment(const Path& path, const std::vector<Alignment>& alignments)
        : path_(path)
        , alignments_(alignments)
    {
        assertValidity();
    }

    uint32_t queryLength() const;
    uint32_t referenceLength() const;
    uint32_t numMatches() const;
    const Path& path() const { return path_; }
    bool overlapsNode(NodeId node_id) const;
    NodeId getNodeIdByIndex(int32_t node_index) const { return path_.getNodeIdByIndex(node_index); }
    std::list<int32_t> getIndexesOfNode(NodeId node_id) const;

    // Removes the specified number of reference bases from the beginning of the alignment while softclipping as many
    // query bases as required
    void shrinkStart(int reference_length);

    // Removes the specified number of reference bases from the end of the alignment while softclipping as many query
    // bases as required
    void shrinkEnd(int reference_length);

    const_iterator begin() const { return alignments_.begin(); }
    const_iterator end() const { return alignments_.end(); }
    const Alignment& front() const { return alignments_.front(); }
    const Alignment& back() const { return alignments_.back(); }
    size_type size() const { return alignments_.size(); }
    const Alignment& operator[](size_t index) const { return alignments_[index]; }
    const std::vector<Alignment>& alignments() const { return alignments_; }
    bool operator==(const GraphAlignment& other) const
    {
        return path_ == other.path_ && alignments_ == other.alignments_;
    }
    bool operator<(const GraphAlignment& other) const;
    std::string generateCigar() const;

    friend std::ostream& operator<<(std::ostream& os, const GraphAlignment& graph_alignment);

private:
    void assertValidity() const;
    Path path_;
    std::vector<Alignment> alignments_;
};

std::ostream& operator<<(std::ostream& os, const GraphAlignment& graph_alignment);
}
