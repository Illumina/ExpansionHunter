// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.

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

#include "graphalign/Operation.hh"

#include "gtest/gtest.h"

using namespace graphtools;

TEST(InitializingOperations, TypicalOperations_QueryAndReferenceSpansObtained)
{
    {
        Operation operation("3M");
        EXPECT_EQ(OperationType::kMatch, operation.type());
        EXPECT_EQ(3u, operation.queryLength());
        EXPECT_EQ(3u, operation.referenceLength());
    }
    {
        Operation operation("4X");
        EXPECT_EQ(OperationType::kMismatch, operation.type());
        EXPECT_EQ(4u, operation.queryLength());
        EXPECT_EQ(4u, operation.referenceLength());
    }
    {
        Operation operation("5D");
        EXPECT_EQ(OperationType::kDeletionFromRef, operation.type());
        EXPECT_EQ(0u, operation.queryLength());
        EXPECT_EQ(5u, operation.referenceLength());
    }
    {
        Operation operation("7I");
        EXPECT_EQ(OperationType::kInsertionToRef, operation.type());
        EXPECT_EQ(7u, operation.queryLength());
        EXPECT_EQ(0u, operation.referenceLength());
    }
    {
        Operation operation("10S");
        EXPECT_EQ(OperationType::kSoftclip, operation.type());
        EXPECT_EQ(10u, operation.queryLength());
        EXPECT_EQ(0u, operation.referenceLength());
    }
    {
        Operation operation("7N");
        EXPECT_EQ(OperationType::kMissingBases, operation.type());
        EXPECT_EQ(7u, operation.queryLength());
        EXPECT_EQ(7u, operation.referenceLength());
    }
}

TEST(EncodingOperation, TypicalOperations_CigarStringObtained)
{
    {
        Operation operation(OperationType::kMatch, 3);
        EXPECT_EQ("3M", operation.generateCigar());
    }
    {
        Operation operation(OperationType::kMismatch, 4);
        EXPECT_EQ("4X", operation.generateCigar());
    }
    {
        Operation operation(OperationType::kDeletionFromRef, 5);
        EXPECT_EQ("5D", operation.generateCigar());
    }
    {
        Operation operation(OperationType::kInsertionToRef, 7);
        EXPECT_EQ("7I", operation.generateCigar());
    }
    {
        Operation operation(OperationType::kSoftclip, 10);
        EXPECT_EQ("10S", operation.generateCigar());
    }
    {
        Operation operation(OperationType::kMissingBases, 7);
        EXPECT_EQ("7N", operation.generateCigar());
    }
}
