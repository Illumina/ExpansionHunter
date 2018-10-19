//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
// Concept: Michael Eberle <meberle@illumina.com>
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

#pragma once

#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "common/count_table.h"
#include "genotyping/RepeatGenotype.hh"

class RepeatFindings
{
public:
    RepeatFindings(
        CountTable countsOfFlankingReads, CountTable countsOfSpanningReads,
        boost::optional<RepeatGenotype> optionalGenotype)
        : countsOfFlankingReads_(countsOfFlankingReads)
        , countsOfSpanningReads_(countsOfSpanningReads)
        , optionalGenotype_(optionalGenotype)
    {
    }

    const CountTable& countsOfFlankingReads() const { return countsOfFlankingReads_; }
    const CountTable& countsOfSpanningReads() const { return countsOfSpanningReads_; }
    const boost::optional<RepeatGenotype>& optionalGenotype() const { return optionalGenotype_; }

    bool operator==(const RepeatFindings& other) const
    {
        return countsOfFlankingReads_ == other.countsOfFlankingReads_
            && countsOfSpanningReads_ == other.countsOfSpanningReads_ && optionalGenotype_ == other.optionalGenotype_;
    }

private:
    CountTable countsOfFlankingReads_;
    CountTable countsOfSpanningReads_;
    boost::optional<RepeatGenotype> optionalGenotype_;
};

std::ostream& operator<<(std::ostream& out, const RepeatFindings& repeatFindings);

using RegionFindings = std::map<std::string, RepeatFindings>;
using SampleFindings = std::map<std::string, RegionFindings>;
