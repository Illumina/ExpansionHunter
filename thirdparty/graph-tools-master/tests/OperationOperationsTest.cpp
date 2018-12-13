//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
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
