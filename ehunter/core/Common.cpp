//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
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

#include "core/Common.hh"

#include <regex>

using std::string;

namespace ehunter
{

Sex decodeSampleSex(const std::string& encoding)
{
    if (encoding == "male")
    {
        return Sex::kMale;
    }
    else if (encoding == "female")
    {
        return Sex::kFemale;
    }
    else
    {
        throw std::invalid_argument(encoding + " is invalid sex; must be either male or female");
    }
}

std::ostream& operator<<(std::ostream& out, Sex sex)
{
    switch (sex)
    {
    case Sex::kFemale:
        out << "Female";
        break;
    case Sex::kMale:
        out << "Male";
        break;
    }

    return out;
}

std::ostream& operator<<(std::ostream& out, ReadType readType)
{
    switch (readType)
    {
    case ReadType::kFlanking:
        out << "FLANKING";
        break;
    case ReadType::kRepeat:
        out << "INREPEAT";
        break;
    case ReadType::kSpanning:
        out << "SPANNING";
        break;
    case ReadType::kOther:
        out << "OTHER";
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, AlleleCount alleleCount)
{
    switch (alleleCount)
    {
    case AlleleCount::kOne:
        out << "One";
        break;
    case AlleleCount::kTwo:
        out << "Two";
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, NumericInterval numericInterval)
{
    out << numericInterval.start() << "-" << numericInterval.end();
    return out;
}

bool isURL(const std::string& path)
{
    static const std::regex url_regex(".*?://.*");
    return std::regex_match(path, url_regex);
}

}
