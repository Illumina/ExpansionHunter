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

#include "graphs/linear_mapping.h"

#include "gtest/gtest.h"

using std::string;
using std::vector;

TEST(OperationInitialization, EncodingOfMatchOperation_OperationCreated) {
  Operation operation("3M", "ATC", "ATC");
  Operation expected_operation('M', 3, "ATC", "ATC");
  ASSERT_EQ(expected_operation, operation);
}

TEST(OperationInitialization, TypicalIncorrectEncodings_ExceptionThrown) {
  EXPECT_ANY_THROW(Operation("4M", "AAAA", "ATCG"));
  EXPECT_ANY_THROW(Operation("4M", "AAAA", "ATC"));
  EXPECT_ANY_THROW(Operation("4M", "AAA", "AAA"));

  EXPECT_ANY_THROW(Operation("4N", "NNN", "NNN"));
  EXPECT_ANY_THROW(Operation("3N", "NN", "NNN"));
  EXPECT_ANY_THROW(Operation("2N", "NT", "NT"));

  EXPECT_ANY_THROW(Operation("2X", "AT", "TT"));
  EXPECT_ANY_THROW(Operation("2X", "AT", "A"));

  EXPECT_ANY_THROW(Operation("4D", "AAA", ""));
  EXPECT_ANY_THROW(Operation("4D", "", ""));

  EXPECT_ANY_THROW(Operation("2I", "AA", "T"));

  EXPECT_ANY_THROW(Operation("2S", "TTT", ""));
  EXPECT_ANY_THROW(Operation("2S", "TT", "T"));
}

TEST(GettingOperationSpans, TypicalOperations_QueryAndReferenceSpansObtained) {
  {
    Operation operation("3M", "AAA", "AAA");
    EXPECT_EQ((int32_t)3, operation.QuerySpan());
    EXPECT_EQ((int32_t)3, operation.ReferenceSpan());
  }
  {
    Operation operation("4X", "AAAA", "TTTT");
    EXPECT_EQ((int32_t)4, operation.QuerySpan());
    EXPECT_EQ((int32_t)4, operation.ReferenceSpan());
  }
  {
    Operation operation("5D", "", "AAAAA");
    EXPECT_EQ((int32_t)0, operation.QuerySpan());
    EXPECT_EQ((int32_t)5, operation.ReferenceSpan());
  }
  {
    Operation operation("7I", "AAAAAAA", "");
    EXPECT_EQ((int32_t)7, operation.QuerySpan());
    EXPECT_EQ((int32_t)0, operation.ReferenceSpan());
  }
  {
    Operation operation("10S", "AAAAAAAAAA", "");
    EXPECT_EQ((int32_t)10, operation.QuerySpan());
    EXPECT_EQ((int32_t)0, operation.ReferenceSpan());
  }
  {
    Operation operation("7N", "NNNNNNN", "NNNNNNN");
    EXPECT_EQ((int32_t)7, operation.QuerySpan());
    EXPECT_EQ((int32_t)7, operation.ReferenceSpan());
  }
  {
    Operation operation("3N", "NCN", "CNN");
    EXPECT_EQ((int32_t)3, operation.QuerySpan());
    EXPECT_EQ((int32_t)3, operation.ReferenceSpan());
  }
}

TEST(EncodingOperation, TypicalOperations_CigarStringObtained) {
  {
    Operation operation("3M", "AAA", "AAA");
    EXPECT_EQ("3M", operation.GetCigarString());
  }
  {
    Operation operation("4X", "AAAA", "TTTT");
    EXPECT_EQ("4X", operation.GetCigarString());
  }
  {
    Operation operation("5D", "", "AAAAA");
    EXPECT_EQ("5D", operation.GetCigarString());
  }
  {
    Operation operation("7I", "AAAAAAA", "");
    EXPECT_EQ("7I", operation.GetCigarString());
  }
  {
    Operation operation("10S", "AAAAAAAAAA", "");
    EXPECT_EQ("10S", operation.GetCigarString());
  }
  {
    Operation operation("7N", "NNNNNNN", "NNNNNNN");
    EXPECT_EQ("7N", operation.GetCigarString());
  }
  {
    Operation operation("3N", "NCN", "CNN");
    EXPECT_EQ("3N", operation.GetCigarString());
  }
}

TEST(MappingInitialization, TypicalCigarString_MappingCreated) {
  // query: ---TTCGTT--TTGGGTCCCCCCCCCC
  //           ||| ||  ||   |
  //   ref: CCCTTCCNNAATT---T----------
  string query = "TTCGTTTTGGGTCCCCCCCCCC";
  string reference = "CCCTTCCNNAATTT";

  Mapping mapping(3, "3M1X2N2D2M3I1M10S", query, reference);
  const vector<Operation> operations = {
      Operation('M', 3, "TTC", "TTC"), Operation('X', 1, "G", "C"),
      Operation('N', 2, "TT", "NN"),   Operation('D', 2, "", "AA"),
      Operation('M', 2, "TT", "TT"),   Operation('I', 3, "GGG", ""),
      Operation('M', 1, "T", "T"),     Operation('S', 10, "CCCCCCCCCC", "")};

  Mapping expected_mapping(3, operations);
  ASSERT_EQ(expected_mapping, mapping);
}

TEST(GettingMappingSpans, TypicalMapping_QueryAndReferenceSpansObtained) {
  Mapping mapping(3, "3M1X2M2D2M3I1M10S", "TTCGTTTTGGGTCCCCCCCCCC",
                  "CCCTTCCTTAATTT");
  EXPECT_EQ((int32_t)22, mapping.QuerySpan());
  EXPECT_EQ((int32_t)11, mapping.ReferenceSpan());
}

TEST(GettingMappingSeqs, TypicalMapping_QueryAndReferenceObtained) {
  Mapping mapping(3, "3M1X2M2D2M3I1M10S", "TTCGTTTTGGGTCCCCCCCCCC",
                  "CCCTTCCTTAATTT");
  EXPECT_EQ("TTCGTTTTGGGT", mapping.Query());
  EXPECT_EQ("TTCCTTAATTT", mapping.Reference());
}

TEST(EncodingMapping, TypicalMapping_CigarStringObtained) {
  string query = "TTCGTTTTGGGTCCCCCCCCCC";
  string reference = "CCCTTCCNNAATTT";
  const std::string cigar_string = "3M1X2N2D2M3I1M10S";

  Mapping mapping(3, cigar_string, query, reference);

  EXPECT_EQ(cigar_string, mapping.GetCigarString());
}
