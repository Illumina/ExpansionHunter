//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Nima Mousavi <mousavi@ucsd.edu>,
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "common/GenomicDistanceArray.hh"

#include <stdexcept>

using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

    std::ostream& operator<<(std::ostream& out, const GenomicDistanceArray& gen_dist_array)
    {
        string encoding;

        for (int32_t element : gen_dist_array.getElements())
        {
            if (!encoding.empty())
            {
                encoding += ", ";
            }

            encoding += to_string(element);
        }

        if (encoding.empty())
        {
            encoding = "";
        }

        out << encoding;

        return out;
    }

}
