//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
