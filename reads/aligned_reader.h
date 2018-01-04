//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>
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

#include "reads/read.h"

namespace reads {

/**
 * @brief Provides access to CRAM/BAM files
 *
 */
class AlignedReader {
 public:
  AlignedReader(const std::string &bam_path, const std::string &reference);
  AlignedReader(AlignedReader &&) noexcept;
  AlignedReader &operator=(AlignedReader &&) noexcept;
  AlignedReader(const AlignedReader &) = delete;
  AlignedReader &operator=(const AlignedReader &) = delete;

  ~AlignedReader();

 private:
  struct Impl;
  std::unique_ptr<Impl> pimpl_;
};

}  // namespace reads