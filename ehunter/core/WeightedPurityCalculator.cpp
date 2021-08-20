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

#include "core/WeightedPurityCalculator.hh"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <unordered_map>

#include <boost/algorithm/string.hpp>

#include "graphutils/SequenceOperations.hh"

using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{
namespace irrdetection
{
using BaseCode = unsigned char;
// const int kMaxBaseCode = 15;

// Core base codes
const BaseCode A = 0;
const BaseCode a = 1;
const BaseCode C = 2;
const BaseCode c = 3;
const BaseCode G = 4;
const BaseCode g = 5;
const BaseCode T = 6;
const BaseCode t = 7;
const BaseCode X = 8;

// Degenerate base codes
const BaseCode B = 9;
const BaseCode D = 10;
const BaseCode H = 11;
const BaseCode K = 12;
const BaseCode M = 13;
const BaseCode N = 14;
const BaseCode R = 15;
const BaseCode S = 16;
const BaseCode V = 17;
const BaseCode W = 18;
const BaseCode Y = 19;

const int kMaxQueryBaseCode = 8;
const int kMaxReferenceBaseCode = 19;

const int maxBaseAscii = 255;

const std::array<BaseCode, maxBaseAscii + 1> kReferenceBaseEncodingTable
    = { X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, A, B, C, D, X, X, G, H, X,
        X, K, X, M, N, X, X, X, R, S, T, X, V, W, X, Y, X, X, X, X, X, X, X, A, X, C, X, X, X, G, X, X, X, X, X, X, X,
        X, X, X, X, X, T, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X };

const std::array<BaseCode, maxBaseAscii + 1> kQueryBaseEncodingTable
    = { X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, A, X, C, X, X, X, G, X, X,
        X, X, X, X, X, X, X, X, X, X, T, X, X, X, X, X, X, X, X, X, X, X, X, a, X, c, X, X, X, g, X, X, X, X, X, X, X,
        X, X, X, X, X, t, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X };

const std::array<std::array<double, kMaxQueryBaseCode + 1>, kMaxReferenceBaseCode + 1>
    kReferenceQueryCodeScoreLookupTable = {
        // clang-format off
          //    A    a     C    c     G    g     T    t     X
             {{ 1.0, 1.0, -1.0, 0.5, -1.0, 0.5, -1.0, 0.5, -1.0}, // A
              { 1.0, 1.0, -1.0, 0.5, -1.0, 0.5, -1.0, 0.5, -1.0}, // a
              {-1.0, 0.5,  1.0, 1.0, -1.0, 0.5, -1.0, 0.5, -1.0}, // C
              {-1.0, 0.5,  1.0, 1.0, -1.0, 0.5, -1.0, 0.5, -1.0}, // c
              {-1.0, 0.5, -1.0, 0.5,  1.0, 1.0, -1.0, 0.5, -1.0}, // G
              {-1.0, 0.5, -1.0, 0.5,  1.0, 1.0, -1.0, 0.5, -1.0}, // g
              {-1.0, 0.5, -1.0, 0.5, -1.0, 0.5,  1.0, 1.0, -1.0}, // T
              {-1.0, 0.5, -1.0, 0.5, -1.0, 0.5,  1.0, 1.0, -1.0}, // t
              {-1.0, 0.5, -1.0, 0.5, -1.0, 0.5, -1.0, 0.5, -1.0}, // X
              {-1.0, 0.5,  1.0, 1.0,  1.0, 1.0,  1.0, 1.0, -1.0}, // B
              { 1.0, 1.0, -1.0, 0.5,  1.0, 1.0,  1.0, 1.0, -1.0}, // D
              { 1.0, 1.0,  1.0, 1.0, -1.0, 0.5,  1.0, 1.0, -1.0}, // H
              {-1.0, 0.5, -1.0, 0.5,  1.0, 1.0,  1.0, 1.0, -1.0}, // K
              { 1.0, 1.0,  1.0, 1.0, -1.0, 0.5, -1.0, 0.5, -1.0}, // M
              { 1.0, 1.0,  1.0, 1.0,  1.0, 1.0,  1.0, 1.0, -1.0}, // N
              { 1.0, 1.0, -1.0, 0.5,  1.0, 1.0, -1.0, 0.5, -1.0}, // R
              {-1.0, 0.5,  1.0, 1.0,  1.0, 1.0, -1.0, 0.5, -1.0}, // S
              { 1.0, 1.0,  1.0, 1.0,  1.0, 1.0, -1.0, 0.5, -1.0}, // V
              { 1.0, 1.0, -1.0, 0.5, -1.0, 0.5,  1.0, 1.0, -1.0}, // W
              {-1.0, 0.5,  1.0, 1.0, -1.0, 0.5,  1.0, 1.0, -1.0}} // Y
        // clang-format on
    };

static inline double scoreBases(char referenceBase, char queryBase)
{
    return kReferenceQueryCodeScoreLookupTable[kReferenceBaseEncodingTable[referenceBase]]
                                              [kQueryBaseEncodingTable[queryBase]];
}
}

WeightedPurityCalculator::WeightedPurityCalculator(const std::string& repeatUnit)
{
    repeatUnits_ = computeCircularPermutations(repeatUnit);
    const string repeatUnitRc = graphtools::reverseComplement(repeatUnit);
    auto permutationsRc = computeCircularPermutations(repeatUnitRc);
    repeatUnits_.insert(repeatUnits_.end(), permutationsRc.begin(), permutationsRc.end());
}

double WeightedPurityCalculator::score(const string& querySequence) const
{
    vector<double> scores;
    for (const auto& repeatUnit : repeatUnits_)
    {
        scores.push_back(score(repeatUnit, querySequence));
    }

    const double weightedPurity
        = *std::max_element(scores.begin(), scores.end()) / static_cast<double>(querySequence.length());
    return weightedPurity;
}

double WeightedPurityCalculator::score(const std::string& repeatUnit, const std::string& querySequence) const
{
    double score = 0;
    for (int position = 0; position != static_cast<int>(querySequence.length()); ++position)
    {
        const int offset = position % repeatUnit.length();
        score += irrdetection::scoreBases(repeatUnit[offset], querySequence[position]);
    }

    return score;
}

vector<string> WeightedPurityCalculator::computeCircularPermutations(string sequence) const
{
    vector<string> permutations = { sequence };
    for (int rotationNum = 0; rotationNum != static_cast<int>(sequence.length()) - 1; ++rotationNum)
    {
        std::rotate(sequence.begin(), sequence.begin() + 1, sequence.end());
        permutations.push_back(sequence);
    }

    return permutations;
}

}
