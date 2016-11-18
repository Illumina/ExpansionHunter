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

#include "include/genomic_region.h"

#include <boost/algorithm/string/split.hpp>
using boost::algorithm::split;
#include <boost/algorithm/string/classification.hpp>
using boost::algorithm::is_any_of;
#include <boost/lexical_cast.hpp>
using boost::lexical_cast;

#include <vector>
using std::vector;
#include <string>
using std::string;
#include <iostream>
using std::endl;
using std::cerr;
#include <sstream>
using std::ostream;
using std::istream;
#include <stdexcept>

/*****************************************************************************/

Region::Region() : chrom_("chr0"), start_(0), end_(0) {}

/*****************************************************************************/

Region::Region(const string& chrom, size_t start, size_t end,
               const string& label)
    : chrom_(chrom), start_(start), end_(end), label_(label) {}

/*****************************************************************************/

Region::Region(const string& encoding, const string& label) : label_(label) {
  vector<string> components;
  split(components, encoding, is_any_of(":-"));

  if (components.size() != 3) {
    throw std::logic_error("Unexpected range format: " + encoding);
  }

  chrom_ = components[0];
  start_ = lexical_cast<size_t>(components[1]);
  end_ = lexical_cast<size_t>(components[2]);
}

/*****************************************************************************/

bool Region::operator<(const Region& other_region) const {
  if (chrom_ != other_region.chrom_) {
    return chrom_ < other_region.chrom_;
  }

  if (start_ != other_region.start_) {
    return start_ < other_region.start_;
  }

  return end_ < other_region.end_;
}

/*****************************************************************************/

bool Region::Overlaps(const Region& other_region) const {
  if (chrom_ != other_region.chrom_) {
    return false;
  }

  const size_t left_bound =
      start_ > other_region.start_ ? start_ : other_region.start_;
  const size_t right_bound =
      end_ < other_region.end_ ? end_ : other_region.end_;

  return left_bound <= right_bound;
}

/*****************************************************************************/

// Returns the range extended by flankSize upstream and downstream.
// NOTE: The right boundary of the extended region may stick past chromosome
// end.
Region Region::Extend(size_t extension_len) const {
  const size_t new_start =
      start_ > extension_len ? (start_ - extension_len) : 1;
  const size_t new_end = end_ + extension_len;
  return Region(chrom_, new_start, new_end);
}

/*****************************************************************************/

const string Region::AsString() const {
  std::ostringstream ostrm;
  ostrm << *this;
  return ostrm.str();
}

/*****************************************************************************/

istream& operator>>(istream& istrm, Region& region) {
  string encoding;
  istrm >> encoding;
  region = Region(encoding);

  return istrm;
}

/*****************************************************************************/

ostream& operator<<(ostream& ostrm, const Region& region) {
  ostrm << region.chrom_ << ':' << region.start_;

  if (region.end_ != region.start_) {
    ostrm << '-' << region.end_;
  }

  if (!region.label_.empty()) {
    ostrm << " " << region.label_;
  }

  return ostrm;
}
