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

#include "reads/read.h"

#include "gtest/gtest.h"

using namespace reads;

TEST(ReadInitialization, TypicalCoreInfo_CoreInfoAddedToRead) {
  Read read;
  read.SetCoreInfo("frag1", "ATTC", "????");
  EXPECT_EQ("frag1", read.FragmentId());
  EXPECT_EQ("ATTC", read.Bases());
  EXPECT_EQ("????", read.Quals());
}

TEST(ReadInitialization, UnsetCoreInfo_ExceptionThrownOnAccess) {
  Read read;
  EXPECT_ANY_THROW(read.FragmentId());
  EXPECT_ANY_THROW(read.Bases());
  EXPECT_ANY_THROW(read.Quals());
}

TEST(ReadInitialization, BasesAndQualsOfUnequalLength_ExceptionThrown) {
  Read read;
  EXPECT_ANY_THROW(read.SetCoreInfo("frag1", "ATT", "?"));
}

TEST(ReadInitialization, TypicalGraphMapping_GraphMappingAddedToRead) {
  Read read;
  read.SetCoreInfo("frag1", "ATTC", "????");
  GraphMappingPtr graph_mapping_ptr(new GraphMapping());
  read.SetCanonicalMapping(std::move(graph_mapping_ptr));
}

TEST(ReadInitialization, UnsertCanonicalMapping_ExceptionThrownOnAccess) {
  Read read;
  read.SetCoreInfo("frag1", "ATTC", "????");
  ASSERT_ANY_THROW(read.CanonicalMapping());
}
