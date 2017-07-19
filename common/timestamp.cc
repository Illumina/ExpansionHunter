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

#include "common/timestamp.h"

std::string TimeStamp() {
  std::time_t now = time(0);
  const size_t timestamp_size = 80;
  char timestamp_buf[timestamp_size];
  std::strftime(timestamp_buf, timestamp_size, "%FT%T", std::localtime(&now));
  return std::string(timestamp_buf);
}