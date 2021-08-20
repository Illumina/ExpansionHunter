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
