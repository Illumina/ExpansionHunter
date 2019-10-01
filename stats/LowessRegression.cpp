//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Xiao Chen <xchen2@illumina.com>
//         Egor Dolzhenko <edolzhenko@illumina.com>
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
/*
 * The implementation closely follows <http://netlib.org/go/lowess.f>,
 * which is a part of the netlib library distributed under license
 * <http://www.netlib.org/math/license.html>
 * original author:
 * wsc@research.bell-labs.com Mon Dec 30 16:55 EST 1985
 * W. S. Cleveland
 * Bell Laboratories
 * Murray Hill NJ 07974
 */

#include "stats/LowessRegression.hh"
#include <algorithm>
#include <cassert>
#include <cmath>
using std::vector;

namespace ehunter
{

void LowessRegression::regression(
    const std::vector<double>& inputXValues, const std::vector<double>& inputYValues,
    std::vector<double>& fittedYValues, std::vector<double>& robustnessWeights, std::vector<double>& residuals) const
{
    int vectorSize = inputXValues.size();
    assert(vectorSize == (int)inputYValues.size());
    std::vector<double> weights(vectorSize);

    if (vectorSize == 1)
    {
        fittedYValues[0] = inputYValues[0];
        return;
    }

    // Number of points bounded by 2 and the vector size
    int NumberOfPointsInNeighborhood = std::max(2, std::min(vectorSize, (int)(smoothingSpan_ * vectorSize)));

    for (int iteration = 0; iteration != numberIteration_ + 1; iteration++)
    {
        bool fitIsOkay;
        int leftBoundary = 0;
        int rightBoundary = NumberOfPointsInNeighborhood - 1;
        int previouslyEstimatedPointIndex = -1;
        int currentPointIndex = 0;
        while (previouslyEstimatedPointIndex < vectorSize - 1)
        {
            // 1. determine neighborhood for fitting
            findNeighbor(inputXValues, currentPointIndex, leftBoundary, rightBoundary);
            // 2. fit
            bool notFirstIteration = (iteration > 0);
            double currentXValue = inputXValues[currentPointIndex];
            double fittedYValue;

            localFit(
                inputXValues, inputYValues, currentXValue, fittedYValue, leftBoundary, rightBoundary, weights,
                notFirstIteration, robustnessWeights, fitIsOkay);

            if (!fitIsOkay)
            {
                fittedYValues[currentPointIndex] = inputYValues[currentPointIndex];
            }
            else
            {
                fittedYValues[currentPointIndex] = fittedYValue;
            }
            // 3. interpolation for points skipped by delta
            if (previouslyEstimatedPointIndex < currentPointIndex - 1)
            {
                interpolate(inputXValues, fittedYValues, currentPointIndex, previouslyEstimatedPointIndex);
            }
            // 4. update indices
            updateIndices(inputXValues, fittedYValues, currentPointIndex, previouslyEstimatedPointIndex);
        }
        // 5. calculate residuals and update weights for the next iteration
        for (int posIndex = 0; posIndex != vectorSize; posIndex++)
        {
            residuals[posIndex] = inputYValues[posIndex] - fittedYValues[posIndex];
        }
        // compute robustness weights except during the last iteration
        if (iteration == numberIteration_)
        {
            return;
        }
        updateWeights(robustnessWeights, residuals, vectorSize);
    }
}

void LowessRegression::findNeighbor(
    const std::vector<double>& inputXValues, int currentPointIndex, int& leftBoundary, int& rightBoundary) const
{
    while (rightBoundary < (int)inputXValues.size() - 1)
    {
        double distanceToLeftBoundary = inputXValues[currentPointIndex] - inputXValues[leftBoundary];
        double distanceToRightBoundary = inputXValues[rightBoundary + 1] - inputXValues[currentPointIndex];
        if (distanceToLeftBoundary <= distanceToRightBoundary)
        {
            return;
        }
        leftBoundary++;
        rightBoundary++;
    }
}

// the original lowest function
void LowessRegression::localFit(
    const std::vector<double>& inputXValues, const std::vector<double>& inputYValues, double currentXValue,
    double& fittedYValue, int leftBoundary, int rightBoundary, std::vector<double>& weights, bool notFirstIteration,
    std::vector<double>& robustnessWeights, bool& fitIsOkay) const
{
    double maxDistance
        = std::max(currentXValue - inputXValues[leftBoundary], inputXValues[rightBoundary] - currentXValue);

    int rightmostPointIndex = 0;
    calculateWeights(
        inputXValues, maxDistance, currentXValue, leftBoundary, weights, notFirstIteration, robustnessWeights,
        fitIsOkay, rightmostPointIndex);

    if (fitIsOkay)
    {
        fitBasedOnWeights(
            inputXValues, inputYValues, fittedYValue, currentXValue, leftBoundary, rightmostPointIndex, weights,
            maxDistance);
    }
}

void LowessRegression::interpolate(
    const std::vector<double>& inputXValues, std::vector<double>& fittedYValues, int currentPointIndex,
    int previouslyEstimatedPointIndex) const
{
    double ratioDenominator = inputXValues[currentPointIndex] - inputXValues[previouslyEstimatedPointIndex];
    for (int posIndex = previouslyEstimatedPointIndex + 1; posIndex != currentPointIndex; posIndex++)
    {
        double ratio = (inputXValues[posIndex] - inputXValues[previouslyEstimatedPointIndex]) / ratioDenominator;
        fittedYValues[posIndex]
            = ratio * fittedYValues[currentPointIndex] + (1 - ratio) * fittedYValues[previouslyEstimatedPointIndex];
    }
}

void LowessRegression::updateIndices(
    const std::vector<double>& inputXValues, std::vector<double>& fittedYValues, int& currentPointIndex,
    int& previouslyEstimatedPointIndex) const
{
    previouslyEstimatedPointIndex = currentPointIndex;
    for (currentPointIndex = previouslyEstimatedPointIndex + 1; currentPointIndex != (int)inputXValues.size();
         currentPointIndex++)
    {
        if (inputXValues[currentPointIndex] > inputXValues[previouslyEstimatedPointIndex] + deltaSkippingParameter_)
        {
            break;
        }
        if (inputXValues[currentPointIndex] == inputXValues[previouslyEstimatedPointIndex])
        {
            fittedYValues[currentPointIndex] = fittedYValues[previouslyEstimatedPointIndex];
            previouslyEstimatedPointIndex = currentPointIndex;
        }
    }
    currentPointIndex = std::max(previouslyEstimatedPointIndex + 1, currentPointIndex - 1);
}

void LowessRegression::updateWeights(
    std::vector<double>& robustnessWeights, const std::vector<double>& residuals, int vectorSize) const
{
    std::vector<double> absoluteResiduals;
    for (int posIndex = 0; posIndex != vectorSize; posIndex++)
    {
        absoluteResiduals.emplace_back(std::abs(residuals[posIndex]));
    }

    std::sort(absoluteResiduals.begin(), absoluteResiduals.end());
    int medianPosition1 = vectorSize / 2;
    int medianPosition2 = vectorSize - vectorSize / 2 - 1;
    // 6 times the median of the absolute residuals
    double sixMedians = 3 * (absoluteResiduals[medianPosition1] + absoluteResiduals[medianPosition2]);
    double upperBound = 0.999 * sixMedians;
    double lowerBound = 0.001 * sixMedians;
    for (int posIndex = 0; posIndex != vectorSize; posIndex++)
    {
        double absoluteResidual = std::abs(residuals[posIndex]);
        if (absoluteResidual <= lowerBound)
        {
            robustnessWeights[posIndex] = 1;
        }
        else if (absoluteResidual > upperBound)
        {
            robustnessWeights[posIndex] = 0;
        }
        else
        {
            robustnessWeights[posIndex] = pow((1 - pow((absoluteResidual / sixMedians), 2)), 2);
        }
    }
}

void LowessRegression::calculateWeights(
    const std::vector<double>& inputXValues, double maxDistance, double currentXValue, int leftBoundary,
    std::vector<double>& weights, bool notFirstIteration, std::vector<double>& robustnessWeights, bool& fitIsOkay,
    int& rightmostPointIndex) const
{
    double sumOfWeights = 0;
    int posIndex;
    // compute weights based on distance
    for (posIndex = leftBoundary; posIndex != (int)inputXValues.size(); posIndex++)
    {
        weights[posIndex] = 0;
        double absoluteDistance = std::abs(inputXValues[posIndex] - currentXValue);
        if (absoluteDistance <= 0.999 * maxDistance)
        {
            if (absoluteDistance > 0.001 * maxDistance)
            {
                double distanceRatio = absoluteDistance / maxDistance;
                weights[posIndex] = pow((1 - pow(distanceRatio, 3)), 3);
            }
            else
            {
                weights[posIndex] = 1.0;
            }
            if (notFirstIteration)
            {
                weights[posIndex] *= robustnessWeights[posIndex];
            }
            sumOfWeights += weights[posIndex];
        }
        // ties can happen here so that now we have passed rightBoundary
        else if (inputXValues[posIndex] > currentXValue)
        {
            break;
        }
    }
    // rightmost point can be greater than rightBoundary because of ties
    rightmostPointIndex = posIndex - 1;

    if (sumOfWeights <= 0)
    {
        fitIsOkay = false;
    }
    else
    {
        fitIsOkay = true;
        for (posIndex = leftBoundary; posIndex != rightmostPointIndex + 1; posIndex++)
        {
            // normalize so that the sum of all weights is 1
            weights[posIndex] /= sumOfWeights;
        }
    }
}

void LowessRegression::fitBasedOnWeights(
    const std::vector<double>& inputXValues, const std::vector<double>& inputYValues, double& fittedYValue,
    double currentXValue, int leftBoundary, int rightmostPointIndex, std::vector<double>& weights,
    double maxDistance) const
{
    if (maxDistance > 0)
    {
        double weightedCenterOfX = 0;

        for (int posIndex = leftBoundary; posIndex != rightmostPointIndex + 1; posIndex++)
        {
            weightedCenterOfX += weights[posIndex] * inputXValues[posIndex];
        }
        double weightedSumOfSquaredDeviations = 0;
        for (int posIndex = leftBoundary; posIndex != rightmostPointIndex + 1; posIndex++)
        {
            weightedSumOfSquaredDeviations += weights[posIndex] * pow((inputXValues[posIndex] - weightedCenterOfX), 2);
        }
        // determine if points are spread out enough to compute slope
        if (pow(weightedSumOfSquaredDeviations, 0.5)
            > 0.001 * (inputXValues[(int)inputXValues.size() - 1] - inputXValues[0]))
        {
            double ratio = (currentXValue - weightedCenterOfX) / weightedSumOfSquaredDeviations;
            for (int posIndex = leftBoundary; posIndex != rightmostPointIndex + 1; posIndex++)
            {
                weights[posIndex] *= 1 + ratio * (inputXValues[posIndex] - weightedCenterOfX);
            }
        }
    }

    fittedYValue = 0;
    for (int posIndex = leftBoundary; posIndex != rightmostPointIndex + 1; posIndex++)
    {
        fittedYValue += weights[posIndex] * inputYValues[posIndex];
    }
}
}
