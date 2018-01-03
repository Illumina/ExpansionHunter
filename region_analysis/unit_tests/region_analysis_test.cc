//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "usecases/reads_operations.h"

#include "gmock/gmock.h"

#include "reads/read.h"
#include "reads/read_pairs.h"
#include "reads/read_reader.h"

using namespace reads;
using ::testing::DoAll;
using ::testing::Return;
using ::testing::SetArgReferee;
using ::testing::_;

class MockReader : public reads::ReadReader {
 public:
  MOCK_METHOD1(SetRegion, void(const Region &region));
  MOCK_METHOD0(GetRead, reads::ReadPtr(void));
};

TEST(ReadExtraction, TypicalRegion_ReadsExtracted) {
  Region target_region("chr1:1-100");
  MockReader mock_reader;
  ReadPtr read1 = std::make_shared<Read>("frag1", "ATCG", "####");
  read1->SetIsFirstMate(true);
  ReadPtr read2 = std::make_shared<Read>("frag1", "GCTA", "####");
  read2->SetIsFirstMate(false);

  EXPECT_CALL(mock_reader, SetRegion(target_region));
  EXPECT_CALL(mock_reader, GetRead())
      .WillOnce(Return(read1))
      .WillOnce(Return(read2))
      .WillOnce(Return(nullptr));

  ReadPairs read_pairs;
  ExtractReads(target_region, mock_reader, read_pairs);

  ReadPairs expected_read_pairs;
  expected_read_pairs.Add(read1);
  expected_read_pairs.Add(read2);

  ASSERT_EQ(expected_read_pairs, read_pairs);
}