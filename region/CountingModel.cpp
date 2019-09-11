//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
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
//

#include "region/CountingModel.hh"

namespace ehunter
{

void CountingModel::analyze(Read read, boost::optional<Read> mate)
{
    readLengthAccumulator_(read.sequence().length());
    if (mate)
    {
        readLengthAccumulator_(mate->sequence().length());
    }
}

int CountingModel::readCount() const { return boost::accumulators::count(readLengthAccumulator_); }

int CountingModel::meanReadLength() const
{
    if (readCount() == 0)
    {
        return 0;
    }

    return boost::accumulators::mean(readLengthAccumulator_);
}

double CountingModel::depth() const
{
    /*
    const int numberOfStartPositions = leftFlankLength_ + rightFlankLength_ - meanReadLength;
    const double depth = meanReadLength * (static_cast<double>(readCount) / numberOfStartPositions);
 */

    return 0;
}

}
