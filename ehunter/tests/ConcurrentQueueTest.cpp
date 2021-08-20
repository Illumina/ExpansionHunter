//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Chris Saunders <csaunders@illumina.com>
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
//

#include "core/ConcurrentQueue.hh"

#include "gtest/gtest.h"

using namespace ehunter;

TEST(ConcurrentQueueTest, SerialOps)
{
    // Sanity-check non-concurrent operation
    ConcurrentQueue<int> cq;
    cq.push(1);
    cq.push(2);
    int r;
    cq.pop(r);
    EXPECT_FALSE(cq.empty());
    EXPECT_EQ(r, 1);
    cq.pop(r);
    EXPECT_TRUE(cq.empty());
    EXPECT_EQ(r, 2);
}
