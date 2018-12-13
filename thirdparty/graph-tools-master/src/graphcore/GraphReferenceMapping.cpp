//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>
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

#include "graphcore/GraphReferenceMapping.hh"

#include <boost/algorithm/string.hpp>
#include <unordered_map>

#include "graphcore/GraphCoordinates.hh"
#include "graphcore/PathOperations.hh"
#include "graphutils/IntervalList.hh"

using std::string;

namespace graphtools
{

ReferenceInterval::ReferenceInterval(ContigId const contig, Position const start, Position const end)
    : contig(contig)
    , start(start)
    , end(end)
{
    if ((start < 0) || (start > end))
    {
        throw std::logic_error("Invalid Interval");
    }
}
ReferenceInterval ReferenceInterval::makePosition(ContigId const& contig, Position const pos)
{
    return ReferenceInterval(contig, pos, pos);
}

ReferenceInterval ReferenceInterval::parseRegion(std::string const& regionString)
{
    std::vector<std::string> spl;
    boost::split(spl, regionString, boost::is_any_of(":"));
    if (spl.size() != 2)
    {
        throw std::runtime_error("Invalid region string: " + regionString);
    }
    string const contig = spl[0];
    boost::replace_all(spl[1], ",", "");
    boost::split(spl, spl[1], boost::is_any_of("-"));
    if (spl.size() != 2)
    {
        throw std::runtime_error("Invalid region string: " + regionString);
    }
    int const start = std::stoll(spl[0]);
    int const end = std::stoll(spl[1]);
    return ReferenceInterval(contig, start, end);
}

bool operator<(ReferenceInterval const& lhs, ReferenceInterval const& rhs)
{
    return lhs.contig < rhs.contig || (lhs.contig == rhs.contig && lhs.start < rhs.start)
        || (lhs.contig == rhs.contig && lhs.start == rhs.start && lhs.end < rhs.end);
}
bool operator==(ReferenceInterval const& lhs, ReferenceInterval const& rhs)
{
    return lhs.contig == rhs.contig && lhs.start == rhs.start && lhs.end == rhs.end;
}
std::ostream& operator<<(std::ostream& output, ReferenceInterval const& ri)
{
    output << ri.contig << ":" << ri.start << "-" << ri.end;
    return output;
}

int32_t ReferenceInterval::length() const { return end - start; }

NodeReferenceMapping::NodeReferenceMapping(Graph const& graph, NodeId const node, ReferenceInterval const& ref)
    : nodeLength_(graph.nodeSeq(node).length())
    , ref_(ref)
{
    if (nodeLength_ != ref.length())
    {
        throw std::logic_error("Length of node sequence does not match reference map length " + graph.nodeName(node));
    }
}

ReferenceInterval NodeReferenceMapping::map(int32_t const offset) const
{
    if (offset >= nodeLength_)
    {
        throw std::logic_error("Cannot map position outside node");
    }
    return ReferenceInterval::makePosition(ref_.contig, ref_.start + offset);
}

GraphReferenceMapping::GraphReferenceMapping(Graph const* graph)
    : graph_(graph)
{
}

void GraphReferenceMapping::addMapping(NodeId node, ReferenceInterval const& ref)
{
    NodeReferenceMapping const nodeMapping(*graph_, node, ref);
    mappings_.insert({ node, nodeMapping });
}

boost::optional<ReferenceInterval> GraphReferenceMapping::map(NodeId node, int32_t offset) const
{
    if (node >= graph_->numNodes())
    {
        throw std::logic_error("Invalid node " + std::to_string(node));
    }
    auto const& mapping = mappings_.find(node);
    if (mapping != mappings_.end())
    {
        return mapping->second.map(offset);
    }
    return boost::optional<ReferenceInterval>();
}

boost::optional<ReferenceInterval> GraphReferenceMapping::map(Path const& path) const
{
    for (auto const& nodePath : generateSubpathForEachNode(path))
    {
        auto const mapped = map(*nodePath.nodeIds().begin(), nodePath.startPosition());
        if (mapped)
        {
            return mapped;
        }
    }
    return boost::optional<ReferenceInterval>();
}
}
