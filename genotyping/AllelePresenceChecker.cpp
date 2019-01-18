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

#include "genotyping/AllelePresenceChecker.hh"

#include <boost/math/special_functions/gamma.hpp>
#include <numeric>

namespace ehunter
{

static double poissonLogLikelihood(double lambda, double count)
{
    return count * log(lambda) - lambda - boost::math::lgamma(count + 1);
}

AllelePresenceStatus
AllelePresenceChecker::check(double haplotypeDepth, int targetAlleleCount, int otherAlleleCount) const
{
    if (haplotypeDepth <= 0)
    {
        throw std::runtime_error("Haplotype depth must be positive");
    }

    if (targetAlleleCount < 0 || otherAlleleCount < 0)
    {
        throw std::runtime_error("Negative read counts are not allowed");
    }

    const int totalReadCount = targetAlleleCount + otherAlleleCount;
    double ll0 = (totalReadCount > 0) ? poissonLogLikelihood(errorRate_ * totalReadCount, targetAlleleCount) : 0;
    double ll1 = poissonLogLikelihood(haplotypeDepth, targetAlleleCount);
    if (abs(ll0 - ll1) < log(llrThreshold_))
    {
        return AllelePresenceStatus::kUncertain;
    }

    return (ll1 > ll0) ? AllelePresenceStatus::kPresent : AllelePresenceStatus::kAbsent;
}

std::ostream& operator<<(std::ostream& out, AllelePresenceStatus status)
{
    switch (status)
    {
    case AllelePresenceStatus::kAbsent:
        out << "Absent";
        break;
    case AllelePresenceStatus::kPresent:
        out << "Present";
        break;
    case AllelePresenceStatus::kUncertain:
        out << "Uncertain";
        break;
    }

    return out;
}

}
