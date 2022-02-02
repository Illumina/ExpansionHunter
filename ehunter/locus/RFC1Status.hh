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

#pragma once

#include <string>

namespace ehunter
{

enum class RFC1CallType
{
    normal,
    potential_carrier,
    carrier,
    affected
};

inline const char* label(RFC1CallType t)
{
    switch (t)
    {
    case RFC1CallType::normal:
        return "normal";
    case RFC1CallType::potential_carrier:
        return "potential carrier";
    case RFC1CallType::carrier:
        return "carrier";
    case RFC1CallType::affected:
        return "affected";
    default:
        return "unknown";
    }
}

struct RFC1Status
{
    RFC1CallType call;
    std::string description;
};

}
