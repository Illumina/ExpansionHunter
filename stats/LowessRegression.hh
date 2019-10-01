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

#pragma once

#include "common/Common.hh"
#include <vector>

namespace ehunter
{

class LowessRegression
{
public:
    LowessRegression(double smoothingSpan, double deltaSkippingParameter, int numberIteration)
        : smoothingSpan_(smoothingSpan)
        , deltaSkippingParameter_(deltaSkippingParameter)
        , numberIteration_(numberIteration)
    {
    }

    void regression(
        const std::vector<double>& inputXValues, const std::vector<double>& inputYValues,
        std::vector<double>& fittedYValues, std::vector<double>& robustnessWeights,
        std::vector<double>& residuals) const;

private:
    // find the neighboring region to perform local fit
    void findNeighbor(
        const std::vector<double>& inputXValues, int currentPointIndex, int& leftBoundary, int& rightBoundary) const;

    // the original lowest function, local fit at each point
    void localFit(
        const std::vector<double>& inputXValues, const std::vector<double>& inputYValues, double currentXValue,
        double& fittedYValue, int leftBoundary, int rightBoundary, std::vector<double>& weights, bool notFirstIteration,
        std::vector<double>& robustnessWeights, bool& fitIsOkay) const;

    // for points skippped by delta, interpolate
    void interpolate(
        const std::vector<double>& inputXValues, std::vector<double>& fittedYValues, int currentPointIndex,
        int previouslyEstimatedPointIndex) const;

    // update indices after skipping some points controlled by delta
    void updateIndices(
        const std::vector<double>& inputXValues, std::vector<double>& fittedYValues, int& currentPointIndex,
        int& previouslyEstimatedPointIndex) const;

    // update weights according to residuals
    void
    updateWeights(std::vector<double>& robustnessWeights, const std::vector<double>& residuals, int vectorSize) const;

    // calculate weights for each point in the neighborhood
    void calculateWeights(
        const std::vector<double>& inputXValues, double maxDistance, double currentXValue, int leftBoundary,
        std::vector<double>& weights, bool notFirstIteration, std::vector<double>& robustnessWeights, bool& fitIsOkay,
        int& rightmostPointIndex) const;

    // fit based on weighted least squares
    void fitBasedOnWeights(
        const std::vector<double>& inputXValues, const std::vector<double>& inputYValues, double& fittedYValue,
        double currentXValue, int leftBoundary, int rightmostPointIndex, std::vector<double>& weights,
        double maxDistance) const;

    // f, specifies the amount of smoothing; f is the fraction of points used to compute each fitted value
    double smoothingSpan_;
    // on the initial fit and on each of the robust fit iterations locally weighted regression fitted values
    // are computed at points which are spaced, roughly, DELTA apart; then the fitted values at the
    // remaining points are computed using linear interpolation.
    double deltaSkippingParameter_;
    // the number of iterations in the robust fit; if numberIteration_ = 0, the nonrobust fit is returned
    double numberIteration_;
};
}
