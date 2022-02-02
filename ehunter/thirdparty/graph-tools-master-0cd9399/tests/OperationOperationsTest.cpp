//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "graphalign/OperationOperations.hh"

#include "gtest/gtest.h"

using namespace graphtools;

TEST(CheckingConsistency, MatchOperation_ConsistencyChecked)
{
    {
        Operation operation("3M");
        EXPECT_TRUE(checkConsistency(operation, "ATC", "ATC"));
    }

    {
        Operation operation("4M");
        EXPECT_TRUE(checkConsistency(operation, "ATBB", "AtcG"));
    }

    {
        Operation operation("4M");
        EXPECT_FALSE(checkConsistency(operation, "AYAA", "AAAA"));
    }

    {
        Operation operation("4M");
        EXPECT_FALSE(checkConsistency(operation, "ATC", "AAAA"));
    }

    {
        Operation operation("4M");
        EXPECT_FALSE(checkConsistency(operation, "AAA", "AAA"));
    }
}

TEST(CheckingConsistency, MismatchOperation_ConsistencyChecked)
{
    {
        Operation operation("2X");
        EXPECT_TRUE(checkConsistency(operation, "TR", "AT"));
    }

    {
        Operation operation("2X");
        EXPECT_FALSE(checkConsistency(operation, "TT", "AT"));
    }

    {
        Operation operation("2X");
        EXPECT_FALSE(checkConsistency(operation, "A", "AT"));
    }

    {
        Operation operation("1X");
        EXPECT_FALSE(checkConsistency(operation, "W", "T"));
    }
}

TEST(CheckingConsistency, InsertionOperation_ConsistencyChecked)
{
    {
        Operation operation("4I");
        EXPECT_TRUE(checkConsistency(operation, "", "ATTG"));
    }

    {
        Operation operation("2I");
        EXPECT_FALSE(checkConsistency(operation, "T", "AA"));
    }
}

TEST(CheckingConsistency, DeletionOperation_ConsistencyChecked)
{

    {
        Operation operation("3D");
        EXPECT_TRUE(checkConsistency(operation, "TRR", ""));
    }

    {
        Operation operation("4D");
        EXPECT_FALSE(checkConsistency(operation, "", "AAA"));
    }

    {
        Operation operation("4D");
        EXPECT_FALSE(checkConsistency(operation, "", ""));
    }
}

TEST(CheckingConsistency, MissingBasesOperation_ConsistencyChecked)
{

    {
        Operation operation("3N");
        EXPECT_TRUE(checkConsistency(operation, "AAN", "NNN"));
    }

    {
        Operation operation("4N");
        EXPECT_FALSE(checkConsistency(operation, "NNN", "NNN"));
    }

    {
        Operation operation("2N");
        EXPECT_FALSE(checkConsistency(operation, "NT", "NT"));
    }

    {
        Operation operation("3N");
        EXPECT_FALSE(checkConsistency(operation, "NNN", "NNA")); // Reference N means degenerate base, not missing base
    }
}

TEST(CheckingConsistency, SoftclipOperation_ConsistencyChecked)
{

    {
        Operation operation("2S");
        EXPECT_TRUE(checkConsistency(operation, "", "AA"));
    }

    {
        Operation operation("2S");
        EXPECT_FALSE(checkConsistency(operation, "", "TTT"));
    }

    {
        Operation operation("2S");
        EXPECT_FALSE(checkConsistency(operation, "T", "TT"));
    }
}

TEST(SplittingOperations, MatchOperation_Split)
{
    Operation operation("3M");
    OperationPair operation_parts = splitByReferenceLength(operation, 1);

    EXPECT_EQ(Operation("1M"), operation_parts.first);
    EXPECT_EQ(Operation("2M"), operation_parts.second);
}

TEST(SplittingOperations, MismatchOperation_Split)
{
    Operation operation("4X");
    OperationPair operation_parts = splitByReferenceLength(operation, 3);

    EXPECT_EQ(Operation("3X"), operation_parts.first);
    EXPECT_EQ(Operation("1X"), operation_parts.second);
}

TEST(SplittingOperations, MissingBaseOperation_ExceptionThrown)
{
    Operation operation("7N");

    OperationPair operation_parts = splitByReferenceLength(operation, 4);

    EXPECT_EQ(Operation("4N"), operation_parts.first);
    EXPECT_EQ(Operation("3N"), operation_parts.second);
}

TEST(SplittingOperations, DeletionOperation_Split)
{
    Operation operation("5D");
    OperationPair operation_parts = splitByReferenceLength(operation, 2);

    EXPECT_EQ(Operation("2D"), operation_parts.first);
    EXPECT_EQ(Operation("3D"), operation_parts.second);
}

TEST(SplittingOperations, InsertionOperation_ExceptionThrown)
{
    Operation operation("7I");

    ASSERT_ANY_THROW(splitByReferenceLength(operation, 2));
}

TEST(SplittingOperations, SoftclipOperation_ExceptionThrown)
{
    Operation operation("10S");

    ASSERT_ANY_THROW(splitByReferenceLength(operation, 2));
}

TEST(SplittingOperations, InvalidReferenceLength_ExceptionThrown)
{

    Operation operation("3M");
    EXPECT_ANY_THROW(splitByReferenceLength(operation, 0));
    EXPECT_ANY_THROW(splitByReferenceLength(operation, 3));
    EXPECT_ANY_THROW(splitByReferenceLength(operation, 4));
}
