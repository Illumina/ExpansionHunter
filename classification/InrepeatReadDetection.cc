//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "classification/InrepeatReadDetection.hh"

#include <algorithm>
#include <stdexcept>
#include <unordered_map>

#include <boost/algorithm/string.hpp>

#include "graphutils/SequenceOperations.hh"

using std::string;
using std::to_string;

double calculatePeriodicityScore(int period, const string& sequence)
{
    if (period <= 0 || static_cast<int>(sequence.length()) / 2 + 1 <= period)
    {
        throw std::logic_error(to_string(period) + " is not a valid period for " + sequence);
    }

    int numMatches = 0;
    for (int position = 0; position != static_cast<int>(sequence.length()) - period; ++position)
    {
        if (sequence[position] == sequence[position + period])
        {
            ++numMatches;
        }
    }

    const int maxMatchesPossible = sequence.length() - period;
    const double matchFrequency = static_cast<double>(numMatches) / maxMatchesPossible;
    return matchFrequency;
}

static char extractConsensusBase(int offset, int period, const string& bases)
{
    std::unordered_map<char, int> charFrequencies;
    for (int index = offset; index < static_cast<int>(bases.length()); index += period)
    {
        ++charFrequencies[bases[index]];
    }

    char consensusChar = '?';
    int maxFrequency = 0;
    for (const auto& charAndFrequency : charFrequencies)
    {
        const char currentChar = charAndFrequency.first;
        const int currentFrequency = charAndFrequency.second;
        if (currentFrequency > maxFrequency)
        {
            maxFrequency = currentFrequency;
            consensusChar = currentChar;
        }
    }

    return consensusChar;
}

string extractConsensusRepeatUnit(int period, const string& bases)
{
    string repeatUnit;
    for (int offset = 0; offset != period; ++offset)
    {
        repeatUnit += extractConsensusBase(offset, period, bases);
    }

    return repeatUnit;
}

string computeSmallestRepeatUnitUnderCircularPermutation(string repeatUnit)
{
    string minimalRepeatUnit = repeatUnit;
    for (int rotationNum = 0; rotationNum != static_cast<int>(repeatUnit.length()) - 1; ++rotationNum)
    {
        std::rotate(repeatUnit.begin(), repeatUnit.begin() + 1, repeatUnit.end());
        if (repeatUnit < minimalRepeatUnit)
        {
            minimalRepeatUnit = repeatUnit;
        }
    }

    return minimalRepeatUnit;
}

string computeCanonicalRepeatUnit(const string& repeatUnit)
{
    string smallestRepeatUnitInCurrentOrientation = computeSmallestRepeatUnitUnderCircularPermutation(repeatUnit);

    const string reverseComplementedUnit = graphtools::reverseComplement(repeatUnit);
    const string smallestReverseComplementedRepeatUnit
        = computeSmallestRepeatUnitUnderCircularPermutation(reverseComplementedUnit);

    if (smallestRepeatUnitInCurrentOrientation < smallestReverseComplementedRepeatUnit)
    {
        return smallestRepeatUnitInCurrentOrientation;
    }
    else
    {
        return smallestReverseComplementedRepeatUnit;
    }
}

bool checkIfSequenceIsRepetitive(const string& repeatUnit, string sequence)
{
    boost::to_upper(sequence);

    const double minPeriodicityScore = 0.75;
    if (calculatePeriodicityScore(repeatUnit.length(), sequence) < minPeriodicityScore)
    {
        return false;
    }

    const string expectedCanonicalUnit = computeCanonicalRepeatUnit(repeatUnit);
    const string canonicalUnit = computeCanonicalRepeatUnit(extractConsensusRepeatUnit(repeatUnit.length(), sequence));

    return expectedCanonicalUnit == canonicalUnit;
}
