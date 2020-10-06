//
// ExpansionHunter
// Copyright (c) 2020 Illumina, Inc.
//
// Author: Konrad Scheffler <kscheffler@illumina.com>
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

// These utility functions were copied from source code of Strelka Small Variant Caller.

#pragma once

#include <boost/math/special_functions/log1p.hpp>

// returns log(1+x), switches to log1p function when abs(x) is small
static double log1p_switch(const double x)
{
    // TODO Justify this switch point. Related discussion:
    // http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
    static const double smallx_thresh(0.01);

    if (std::abs(x) < smallx_thresh)
    {
        return boost::math::log1p(x);
    }
    else
    {
        return std::log(1 + x);
    }
}

// Returns the equivalent of log(exp(x1)+exp(x2))
static double getLogSum(double x1, double x2)
{
    if (x1 < x2)
        std::swap(x1, x2);
    return x1 + log1p_switch(std::exp(x2 - x1));
}