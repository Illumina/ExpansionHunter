// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.

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

#include "graphalign/GappedAligner.hh"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphalign/LinearAlignmentOperations.hh"
#include "graphcore/PathOperations.hh"

using std::list;
using std::make_pair;
using std::string;
using std::to_string;

namespace graphtools
{

/**
 * Remove prefix of a path if it overlaps multiple nodes
 *
 * @param prefixLength: length of the prefix of the path to evaluate
 * @param[in|out] path: any path
 * @return: the length by which the path's beginning was trimmed
 */
static int32_t removePrefixThatOverlapsMultipleNodes(int32_t prefixLength, Path& path)
{
    if (path.numNodes() == 1)
    {
        return 0;
    }

    const int32_t overlapLengthWithFirstNode = path.getNodeOverlapLengthByIndex(0);

    if (overlapLengthWithFirstNode <= prefixLength)
    {
        if ((int32_t)path.length() >= prefixLength)
        {
            path.shrinkStartBy(prefixLength);
            return prefixLength;
        }
    }

    return 0;
}

/**
 * Remove suffix of a path if it overlaps multiple nodes
 *
 * @param suffixLength: length of the suffix of the path to evaluate
 * @param[in|out] path: any path
 * @return: the length by which the path's end was trimmed
 */
static int32_t removeSuffixThatOverlapsMultipleNodes(int32_t suffixLength, Path& path)
{
    if (path.numNodes() == 1)
    {
        return 0;
    }

    const int32_t overlapLengthWithLastNode = path.getNodeOverlapLengthByIndex(path.numNodes() - 1);

    if (overlapLengthWithLastNode <= suffixLength)
    {
        if ((int32_t)path.length() >= suffixLength)
        {
            path.shrinkEndBy(suffixLength);
            return suffixLength;
        }
    }

    return 0;
}

list<GraphAlignment> GappedGraphAligner::align(const string& query) const
{
    const list<string> kmers = extractKmersFromAllPositions(query, kmer_len_);

    size_t kmer_start_on_query = 0;
    for (const string& kmer : kmers)
    {
        // Initiate alignment from a unique kmer.
        if (kmer_index_.numPaths(kmer) == 1)
        {
            Path kmer_path = kmer_index_.getPaths(kmer).front();
            removeSuffixThatOverlapsMultipleNodes(seed_affix_trim_len_, kmer_path);
            const int32_t num_prefix_bases_trimmed
                = removePrefixThatOverlapsMultipleNodes(seed_affix_trim_len_, kmer_path);
            return extendKmerMatchToFullAlignments(kmer_path, query, kmer_start_on_query + num_prefix_bases_trimmed);
        }
        ++kmer_start_on_query;
    }

    return {};
}

list<GraphAlignment> GappedGraphAligner::extendKmerMatchToFullAlignments(
    Path kmer_path, const string& query, size_t kmer_start_on_query) const
{
    assert(kmer_path.length() > 1);

    // Generate prefix extensions
    list<PathAndAlignment> prefix_extensions;
    size_t query_prefix_len = kmer_start_on_query;
    if (query_prefix_len != 0)
    {
        const string query_prefix = query.substr(0, query_prefix_len);
        Path prefix_seed_path = kmer_path;
        prefix_seed_path.shrinkEndBy(kmer_path.length());
        prefix_extensions = extendAlignmentPrefix(prefix_seed_path, query_prefix, query_prefix_len + padding_len_);
    }
    else
    {
        // Because (a) empty alignments are currently disallowed and (b) we don't want to deal with an empty list of
        // prefix_extensions we create a 1bp prefix artificially.
        query_prefix_len = 1;
        Path prefix_path = kmer_path;
        prefix_path.shrinkEndBy(prefix_path.length() - 1);
        prefix_extensions = { make_pair(prefix_path, Alignment(0, "1M")) };
        kmer_path.shrinkStartBy(1);
    }

    // Generate suffix extensions
    list<PathAndAlignment> suffix_extensions;
    size_t query_suffix_len = query.length() - kmer_path.length() - query_prefix_len;
    if (query_suffix_len != 0)
    {
        const string query_suffix = query.substr(query_prefix_len + kmer_path.length(), query_suffix_len);
        Path suffix_seed_path = kmer_path;
        suffix_seed_path.shrinkStartBy(kmer_path.length());
        suffix_extensions = extendAlignmentSuffix(suffix_seed_path, query_suffix, query_suffix_len + padding_len_);
    }
    else
    {
        // Because (a) empty alignments are currently disallowed and (b) we don't want to deal with an empty list of
        // suffix_extensions we create a 1bp suffix artificially.
        Path suffix_path = kmer_path;
        suffix_path.shrinkStartBy(suffix_path.length() - 1);
        suffix_extensions = { make_pair(suffix_path, Alignment(0, "1M")) };
        kmer_path.shrinkEndBy(1);
    }

    // Merge alignments together
    list<PathAndAlignment> top_paths_and_alignments;
    for (PathAndAlignment& prefix_path_and_alignment : prefix_extensions)
    {
        Path& prefix_path = prefix_path_and_alignment.first;
        Path prefix_plus_kmer_path = concatenatePaths(prefix_path, kmer_path);

        Alignment& prefix_alignment = prefix_path_and_alignment.second;
        Alignment kmer_alignment(prefix_alignment.referenceLength(), to_string(kmer_path.length()) + "M");
        Alignment prefix_plus_kmer_alignment = mergeAlignments(prefix_alignment, kmer_alignment);

        for (PathAndAlignment& suffix_path_and_alignment : suffix_extensions)
        {
            Path& suffix_path = suffix_path_and_alignment.first;
            Alignment& suffix_alignment = suffix_path_and_alignment.second;
            Path full_path = concatenatePaths(prefix_plus_kmer_path, suffix_path);

            suffix_alignment.setReferenceStart(prefix_plus_kmer_path.length());
            Alignment full_alignment = mergeAlignments(prefix_plus_kmer_alignment, suffix_alignment);
            top_paths_and_alignments.push_back(make_pair(full_path, full_alignment));
        }
    }

    list<GraphAlignment> top_graph_alignments;
    for (PathAndAlignment& path_and_alignment : top_paths_and_alignments)
    {
        Path& path = path_and_alignment.first;
        Alignment& alignment = path_and_alignment.second;
        top_graph_alignments.push_back(projectAlignmentOntoGraph(alignment, path));
    }

    top_graph_alignments.sort();
    top_graph_alignments.unique();
    return top_graph_alignments;
}

// bool comparePnA(const PathAndAlignment& left, const PathAndAlignment& right)
//{
//    return left.first.encode() < right.first.encode()
//        || ((left.first.encode() == right.first.encode() && left.first.seq() < right.first.seq())
//            || (left.first.seq() == right.first.seq() && (left.second.generateCigar() <
//            right.second.generateCigar())));
//}
//
// bool pathLess(const PathAndAlignment& left, const PathAndAlignment& right)
//{
//    return left.first.encode() < right.first.encode();
//}
//
// bool pathMatch(const PathAndAlignment& left, const PathAndAlignment& right)
//{
//    return left.first.encode() == right.first.encode();
//}
//
// void traceDiff(
//    const std::list<PathAndAlignment>& left, const std::list<PathAndAlignment>& right, int32_t leftScore,
//    int32_t rightScore)
//{
//    std::vector<PathAndAlignment> vLeft(left.begin(), left.end());
//    std::vector<PathAndAlignment> vRight(right.begin(), right.end());
//
//    std::sort(vLeft.begin(), vLeft.end(), comparePnA);
//    std::sort(vRight.begin(), vRight.end(), comparePnA);
//
//    auto leftIt = vLeft.begin(), rightIt = vRight.begin();
//    for (; vLeft.end() != leftIt || vRight.end() != rightIt;)
//    {
//        if (vLeft.end() == leftIt || (vRight.end() != rightIt && pathLess(*rightIt, *leftIt)))
//        {
//            std::cerr << "missing vs " << rightIt->second.referenceStart() << ":" << rightIt->second.generateCigar()
//                      << "(" << rightScore << ") on " << rightIt->first.encode() << ":" << rightIt->first.seq()
//                      << std::endl;
//            ++rightIt;
//        }
//        else if (vRight.end() == rightIt || (vLeft.end() != leftIt && pathLess(*leftIt, *rightIt)))
//        {
//            std::cerr << leftIt->second.referenceStart() << ":" << leftIt->second.generateCigar() << "(" << leftScore
//                      << ") on " << leftIt->first.encode() << ":" << leftIt->first.seq() << " vs missing" <<
//                      std::endl;
//            ++leftIt;
//        }
//        else
//        {
//            if (!pathMatch(*leftIt, *rightIt))
//            {
//                throw std::logic_error("Unexpected path mismatch");
//            }
//            std::cerr << leftIt->second.referenceStart() << ":" << leftIt->second.generateCigar() << "(" << leftScore
//                      << ") vs " << rightIt->second.referenceStart() << ":" << rightIt->second.generateCigar() << "("
//                      << rightScore << ") on " << rightIt->first.encode() << ":" << rightIt->first.seq() << std::endl;
//            ++leftIt;
//            ++rightIt;
//        }
//    }
//}

list<PathAndAlignment>
GappedGraphAligner::extendAlignmentPrefix(const Path& seed_path, const string& query_piece, size_t extension_len) const
{
    assert(seed_path.length() == 0);

    int32_t top_alignment_score = INT32_MIN;
    list<PathAndAlignment> top_paths_and_alignments
        = aligner_.suffixAlign(seed_path, query_piece, extension_len, top_alignment_score);

    //    traceDiff(dag_paths_and_alignments, top_paths_and_alignments, top_dag_score, top_alignment_score);

    for (PathAndAlignment& path_and_alignment : top_paths_and_alignments)
    {
        Path& path = path_and_alignment.first;
        Alignment& alignment = path_and_alignment.second;
        alignment.setReferenceStart(0);

        const int32_t overhang = path.length() - alignment.referenceLength();
        path.shrinkStartBy(overhang);

        if (!checkConsistency(path_and_alignment.second, path.seq(), query_piece))
        {
            // std::cerr << " alignment = " << path_and_alignment.second << " is inconsistent  with";
            // std::cerr << "path = " << path << " having sequence " << path.seq() << std::endl;
            throw std::logic_error("Inconsistent prefix");
        }
    }

    return top_paths_and_alignments;
}

list<PathAndAlignment>
GappedGraphAligner::extendAlignmentSuffix(const Path& seed_path, const string& query_piece, size_t extension_len) const
{
    assert(seed_path.length() == 0);

    int32_t top_alignment_score = INT32_MIN;
    list<PathAndAlignment> top_paths_and_alignments
        = aligner_.prefixAlign(seed_path, query_piece, extension_len, top_alignment_score);

    //    traceDiff(dag_paths_and_alignments, top_paths_and_alignments, top_dag_score, top_alignment_score);

    for (PathAndAlignment& path_and_alignment : top_paths_and_alignments)
    {
        if (!checkConsistency(path_and_alignment.second, path_and_alignment.first.seq(), query_piece))
        {
            throw std::logic_error("Inconsistent suffix");
        }

        Path& path = path_and_alignment.first;
        Alignment& alignment = path_and_alignment.second;

        const int32_t overhang = path.length() - alignment.referenceLength();
        path.shrinkEndBy(overhang);
    }

    return top_paths_and_alignments;
}
}
