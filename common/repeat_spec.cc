//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>

#include "common/ref_genome.h"
#include "common/repeat_spec.h"
#include "common/timestamp.h"
#include "purity/purity.h"

using std::string;
using std::cerr;
using std::endl;
using std::vector;
using std::map;

using boost::property_tree::ptree;
using boost::lexical_cast;
using boost::algorithm::join;

typedef boost::tokenizer<boost::char_separator<char>> Tokenizer;

const char RepeatSpec::LeftFlankBase() const {
  if (left_flank.empty()) {
    return '.';
  }

  return left_flank[left_flank.size() - 1];
}

RepeatSpec::RepeatSpec(const string &json_path) {
  std::ifstream istrm(json_path.c_str());

  if (!istrm.is_open()) {
    throw std::logic_error("Failed to open region JSON file " + json_path);
  }

  ptree root_node;
  boost::property_tree::read_json(istrm, root_node);

  repeat_id = root_node.get_child("RepeatId").data();
  const string unit = root_node.get_child("RepeatUnit").data();

  const boost::char_separator<char> slash_separator("/");

  Tokenizer tokenizer(unit, slash_separator);
  units = vector<string>(tokenizer.begin(), tokenizer.end());
  units_shifts = shift_units(units);

  is_common_unit_ = false;
  boost::optional<ptree> commonality_node(
      root_node.get_child_optional("CommonUnit"));
  if (commonality_node) {
    const string is_common_unit_encoding = commonality_node->data();
    if (is_common_unit_encoding == "true") {
      is_common_unit_ = true;
    } else if (is_common_unit_encoding == "false") {
      is_common_unit_ = false;
    } else {
      throw std::runtime_error(
          "ERROR: CommonUnit must be either \"true\" or \"false\".");
    }
  }

  const string region_encoding(root_node.get_child("TargetRegion").data());
  target_region = Region(region_encoding);

  boost::optional<ptree> confusion_node(
      root_node.get_child_optional("OffTargetRegions"));

  if (confusion_node) {
    offtarget_regions.clear();
    for (const ptree::value_type &region_node : *confusion_node) {
      assert(region_node.first.empty()); // array elements have no names
      offtarget_regions.push_back(Region(region_node.second.data()));
    }
  }
}

// Fill out prefix and suffix sequences.
bool LoadFlanks(const string &genome_path, double min_wp,
                RepeatSpec *repeat_spec) {
  RefGenome ref_genome(genome_path);
  // Reference repeat flanks should be at least as long as reads.
  const int kFlankLen = 250;

  const Region &repeat_region = repeat_spec->target_region;

  const int64_t left_flank_begin = repeat_region.start() - kFlankLen;
  const int64_t left_flank_end = repeat_region.start() - 1;
  const int64_t right_flank_begin = repeat_region.end() + 1;
  const int64_t right_flank_end = repeat_region.end() + kFlankLen;

  const string left_flank_coords = repeat_region.chrom() + ":" +
                                   lexical_cast<string>(left_flank_begin) +
                                   "-" + lexical_cast<string>(left_flank_end);
  const string right_flank_coords = repeat_region.chrom() + ":" +
                                    lexical_cast<string>(right_flank_begin) +
                                    "-" + lexical_cast<string>(right_flank_end);

  ref_genome.ExtractSeq(left_flank_coords, &repeat_spec->left_flank);
  ref_genome.ExtractSeq(right_flank_coords, &repeat_spec->right_flank);

  // Output prefix, suffix, repeat, and whole locus.
  const string repeat_coords = repeat_region.chrom() + ":" +
                               lexical_cast<string>(left_flank_end + 1) + "-" +
                               lexical_cast<string>(right_flank_begin - 1);
  ref_genome.ExtractSeq(repeat_coords, &repeat_spec->ref_seq);

  string fake_quals = string('P', repeat_spec->ref_seq.length());
  double ref_repeat_wp =
      MatchRepeat(repeat_spec->units, repeat_spec->ref_seq, fake_quals);
  ref_repeat_wp /= repeat_spec->ref_seq.length();
  if (ref_repeat_wp < min_wp) {
    cerr << "[WARNING: reference sequence of " << repeat_spec->repeat_id
         << " repeat (" << repeat_spec->ref_seq
         << ") has low weighed purity score of "
         << lexical_cast<string>(ref_repeat_wp) << "]";
  }

  return true;
}

bool LoadRepeatSpecs(const string &specs_path, const string &genome_path,
                     double min_wp, map<string, RepeatSpec> *repeat_specs) {
  assert(!specs_path.empty());

  const boost::regex regex_json(".*\\.json$");
  boost::filesystem::path path(specs_path);
  boost::filesystem::directory_iterator end_itr;

  for (boost::filesystem::directory_iterator itr(path); itr != end_itr; ++itr) {
    if (is_regular_file(itr->status())) {
      const string fname = itr->path().filename().string();
      boost::smatch what;

      if (boost::regex_match(fname, what, regex_json)) {
        cerr << TimeStamp() << ",[Loading " << fname << "]" << endl;

        const string json_path = itr->path().string();
        RepeatSpec repeat_spec(json_path);
        LoadFlanks(genome_path, min_wp, &repeat_spec);
        (*repeat_specs)[repeat_spec.repeat_id] = repeat_spec;
      }
    }
  }

  return true;
}
