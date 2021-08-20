//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Roman Petrovski <RPetrovski@illumina.com>
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

#include "graphalign/GappedAligner.hh"

#include <algorithm>
#include <stdexcept>

#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphalign/LinearAlignmentOperations.hh"
#include "graphcore/PathOperations.hh"

using boost::optional;
using std::list;
using std::logic_error;
using std::make_pair;
using std::string;
using std::to_string;

namespace graphtools
{

/**
 * Trim prefix of a path if it is close to node edge
 *
 * @param requested_trim_len: length of the prefix to attempt to trim by
 * @param min_path_len: minimal length of the trimmed path
 * @param[in|out] path: any path
 * @return: actual trim length
 */
static int32_t trimPrefixNearNodeEdge(int32_t requested_trim_len, int32_t min_path_len, Path& path)
{
    if (path.numNodes() == 1 || static_cast<int32_t>(path.length()) <= min_path_len)
    {
        return 0;
    }

    const int32_t overlap_len_with_first_node = path.getNodeOverlapLengthByIndex(0);
    if (overlap_len_with_first_node > requested_trim_len)
    {
        return 0;
    }

    const bool can_be_fully_trimmed = static_cast<int32_t>(path.length()) >= requested_trim_len + min_path_len;
    const int32_t actual_trim_len = can_be_fully_trimmed ? requested_trim_len : path.length() - min_path_len;

    path.shrinkStartBy(actual_trim_len);
    return actual_trim_len;
}

/**
 * Trim suffix of a path if it is close to node edge
 *
 * @param requested_trim_len: length of the suffix to attempt to trim by
 * @param min_path_len: minimal length of the shrank path
 * @param[in|out] path: any path
 * @return: the length by which the path's end was trimmed
 */
static int32_t trimSuffixNearNodeEdge(int32_t requested_trim_len, int32_t min_path_len, Path& path)
{
    if (path.numNodes() == 1 || static_cast<int32_t>(path.length()) <= min_path_len)
    {
        return 0;
    }

    const int32_t overlap_len_with_last_node = path.getNodeOverlapLengthByIndex(path.numNodes() - 1);
    if (overlap_len_with_last_node > requested_trim_len)
    {
        return 0;
    }

    const bool can_be_fully_trimmed = static_cast<int32_t>(path.length()) >= requested_trim_len + min_path_len;
    const int32_t actual_trim_len = can_be_fully_trimmed ? requested_trim_len : path.length() - min_path_len;

    path.shrinkEndBy(actual_trim_len);
    return actual_trim_len;
}

list<GraphAlignment> GappedGraphAligner::align(const string& query, AlignerSelector& alignerSelector) const
{
    try
    {
        optional<AlignmentSeed> optional_seed = searchForAlignmentSeed(query);

        if (optional_seed)
        {
            Path& seed_path = optional_seed->path;
            int seed_start_on_query = optional_seed->start_on_query;

            const int kMinPathLength = 2;
            trimSuffixNearNodeEdge(seed_affix_trim_len_, kMinPathLength, seed_path);
            const int trimmed_prefix_len = trimPrefixNearNodeEdge(seed_affix_trim_len_, kMinPathLength, seed_path);
            return extendSeedToFullAlignments(
                seed_path, query, seed_start_on_query + trimmed_prefix_len, alignerSelector);
        }
        else
        {
            return {};
        }
    }
    catch (const std::exception& e)
    {
        throw logic_error("Unable to align " + query + ": " + e.what());
    }
}

optional<GappedGraphAligner::AlignmentSeed> GappedGraphAligner::searchForAlignmentSeed(const string& query) const
{
    string upperQuery = query;
    boost::to_upper(upperQuery);

    optional<GappedGraphAligner::AlignmentSeed> optional_seed;

    bool found_multipath_kmer = false;
    size_t kmer_start_position = 0;
    while (kmer_start_position + kmer_len_ <= upperQuery.length())
    {
        const string kmer = upperQuery.substr(kmer_start_position, static_cast<size_t>(kmer_len_));

        // Initiate seed construction from a unique kmer
        auto num_kmer_paths = kmer_index_.numPaths(kmer);
        if (num_kmer_paths > 1)
        {
            found_multipath_kmer = true;
        }

        if (num_kmer_paths == 1)
        {
            const Path kmer_path = kmer_index_.getPaths(kmer).front();
            // This call updates kmer_start_position to the start of the extended path
            Path extended_path = extendPathMatching(kmer_path, upperQuery, kmer_start_position);

            if (!optional_seed || extended_path.length() > optional_seed->path.length())
            {
                optional_seed = AlignmentSeed(extended_path, kmer_start_position);
            }

            kmer_start_position += extended_path.length();
        }
        else
        {
            ++kmer_start_position;
        }
    }

    if (optional_seed || !found_multipath_kmer)
    {
        return optional_seed;
    }

    // If the search for unique kmer failed, consider kmers that correspond to multiple paths
    const int kMaxPathCount = 10;
    kmer_start_position = 0;
    while (kmer_start_position + kmer_len_ <= upperQuery.length())
    {
        const string kmer = upperQuery.substr(kmer_start_position, static_cast<size_t>(kmer_len_));

        const int numPaths = kmer_index_.numPaths(kmer);
        if (0 < numPaths && numPaths <= kMaxPathCount)
        {
            size_t longest_kmer_path_extension = 0;
            size_t kmer_start_position_for_longest_extension = 0;
            for (const Path& kmer_path : kmer_index_.getPaths(kmer))
            {
                size_t kmer_start_position_for_kmer_path = kmer_start_position;
                Path extended_path = extendPathMatching(kmer_path, upperQuery, kmer_start_position_for_kmer_path);

                if (longest_kmer_path_extension < extended_path.length())
                {
                    longest_kmer_path_extension = extended_path.length();
                    kmer_start_position_for_longest_extension = kmer_start_position_for_kmer_path;
                }

                if (!optional_seed || extended_path.length() > optional_seed->path.length())
                {
                    optional_seed = AlignmentSeed(extended_path, kmer_start_position_for_kmer_path);
                }
            }

            kmer_start_position = kmer_start_position_for_longest_extension + longest_kmer_path_extension;
        }
        else
        {
            ++kmer_start_position;
        }
    }

    return optional_seed;
}

list<GraphAlignment> GappedGraphAligner::extendSeedToFullAlignments(
    Path seed_path, const string& query, size_t seed_start_on_query, AlignerSelector& alignerSelector) const
{
    assert(seed_path.length() > 1);

    // Generate prefix extensions
    list<PathAndAlignment> prefix_extensions;
    size_t query_prefix_len = seed_start_on_query;
    if (query_prefix_len != 0)
    {
        const string query_prefix = query.substr(0, query_prefix_len);
        Path prefix_seed_path = seed_path;
        prefix_seed_path.shrinkEndBy(seed_path.length());
        prefix_extensions
            = extendAlignmentPrefix(prefix_seed_path, query_prefix, query_prefix_len + padding_len_, alignerSelector);
    }
    else
    {
        // Because (a) empty alignments are currently disallowed and (b) we don't want to deal with an empty list of
        // prefix_extensions we create a 1bp prefix artificially.
        query_prefix_len = 1;
        Path prefix_path = seed_path;
        prefix_path.shrinkEndBy(prefix_path.length() - 1);
        prefix_extensions = { make_pair(prefix_path, Alignment(0, "1M")) };
        seed_path.shrinkStartBy(1);
    }

    // Generate suffix extensions
    list<PathAndAlignment> suffix_extensions;
    size_t query_suffix_len = query.length() - seed_path.length() - query_prefix_len;
    if (query_suffix_len != 0)
    {
        const string query_suffix = query.substr(query_prefix_len + seed_path.length(), query_suffix_len);
        Path suffix_seed_path = seed_path;
        suffix_seed_path.shrinkStartBy(seed_path.length());
        suffix_extensions
            = extendAlignmentSuffix(suffix_seed_path, query_suffix, query_suffix_len + padding_len_, alignerSelector);
    }
    else
    {
        // Because (a) empty alignments are currently disallowed and (b) we don't want to deal with an empty list of
        // suffix_extensions we create a 1bp suffix artificially.
        Path suffix_path = seed_path;
        suffix_path.shrinkStartBy(suffix_path.length() - 1);
        suffix_extensions = { make_pair(suffix_path, Alignment(0, "1M")) };
        seed_path.shrinkEndBy(1);
    }

    // Merge alignments together
    list<PathAndAlignment> top_paths_and_alignments;
    for (PathAndAlignment& prefix_path_and_alignment : prefix_extensions)
    {
        Path& prefix_path = prefix_path_and_alignment.first;
        Path prefix_plus_seed_path = concatenatePaths(prefix_path, seed_path);

        Alignment& prefix_alignment = prefix_path_and_alignment.second;
        Alignment kmer_alignment(prefix_alignment.referenceLength(), to_string(seed_path.length()) + "M");
        Alignment prefix_plus_kmer_alignment = mergeAlignments(prefix_alignment, kmer_alignment);

        for (PathAndAlignment& suffix_path_and_alignment : suffix_extensions)
        {
            Path& suffix_path = suffix_path_and_alignment.first;
            Alignment& suffix_alignment = suffix_path_and_alignment.second;
            Path full_path = concatenatePaths(prefix_plus_seed_path, suffix_path);

            suffix_alignment.setReferenceStart(prefix_plus_seed_path.length());
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

list<PathAndAlignment> GappedGraphAligner::extendAlignmentPrefix(
    const Path& seed_path, const string& query_piece, size_t extension_len, AlignerSelector& alignerSelector) const
{
    assert(seed_path.length() == 0);

    int32_t top_alignment_score = INT32_MIN;
    list<PathAndAlignment> top_paths_and_alignments
        = alignerSelector.suffixAlign(seed_path, query_piece, extension_len, top_alignment_score);

    for (PathAndAlignment& path_and_alignment : top_paths_and_alignments)
    {
        Path& path = path_and_alignment.first;
        Alignment& alignment = path_and_alignment.second;
        alignment.setReferenceStart(0);

        const int32_t overhang = path.length() - alignment.referenceLength();
        path.shrinkStartBy(overhang);

        if (!checkConsistency(path_and_alignment.second, path.seq(), query_piece))
        {
            throw std::logic_error("Inconsistent prefix");
        }
    }

    return top_paths_and_alignments;
}

list<PathAndAlignment> GappedGraphAligner::extendAlignmentSuffix(
    const Path& seed_path, const string& query_piece, size_t extension_len, AlignerSelector& alignerSelector) const
{
    assert(seed_path.length() == 0);

    int32_t top_alignment_score = INT32_MIN;
    list<PathAndAlignment> top_paths_and_alignments
        = alignerSelector.prefixAlign(seed_path, query_piece, extension_len, top_alignment_score);

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
