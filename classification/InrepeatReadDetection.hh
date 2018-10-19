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

#pragma once
#include <string>

// Given a sequence s and a period p, the periodicity score is defined to be the fraction of bases satisfying
// s(i) = s(i + p)
double calculatePeriodicityScore(int period, const std::string& sequence);

// Determines the consensus repeat unit for a given period
std::string extractConsensusRepeatUnit(int period, const std::string& bases);

// Determines the smallest repeat unit (in lexicographic order) that can be obtained from the given one by performing
// circular permutations
std::string computeSmallestRepeatUnitUnderCircularPermutation(std::string repeatUnit);

// Computes the canonical representation of a given repeat unit; the canonical representation is defined to be the
// smallest repeat unit that could be obtained from the given one by performing circular permutations and reverse
// complement operations
std::string computeCanonicalRepeatUnit(const std::string& repeatUnit);

// A sequence is assumed to be repetitive with the given repeat unit if (a) it has a high periodicity score and (b)
// its consensus repeat unit matches the one provided
bool checkIfSequenceIsRepetitive(const std::string& repeatUnit, std::string sequence);
