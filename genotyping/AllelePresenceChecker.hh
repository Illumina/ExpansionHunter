//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
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

#include <iostream>

#include "common/Common.hh"

namespace ehunter
{

enum class AllelePresenceStatus
{
    kPresent,
    kAbsent,
    kUncertain
};

/*
 * Genotyper checking for presence (>= 1 allele) of a given 'Key' allele
 */
class AllelePresenceChecker
{
public:
    AllelePresenceChecker(double haplotypeDepth, double errorRate = 0.02, double llrThreshold = 10000)
        : haplotypeDepth_(haplotypeDepth)
        , errorRate_(errorRate)
        , llrThreshold_(llrThreshold)
    {
        if (haplotypeDepth <= 0)
        {
            throw std::runtime_error("Haplotype depth must be positive");
        }
        if (errorRate <= 0 || errorRate >= 1)
        {
            throw std::runtime_error("Error rate must be positive and less than 1");
        }
        if (llrThreshold < 0)
        {
            throw std::runtime_error("Likelihood Ratio threshold must be positive");
        }
    }

    AllelePresenceStatus check(int targetAlleleCount, int otherAlleleCount) const;

private:
    // expected depth for one allele
    double haplotypeDepth_;
    // Rate of 'false' key-allele observations
    double errorRate_;
    // If the likelihood ratio threshold in favor of presence or absence
    // is not at least this strong, return no call.
    double llrThreshold_;
};

std::ostream& operator<<(std::ostream& out, AllelePresenceStatus status);

}
