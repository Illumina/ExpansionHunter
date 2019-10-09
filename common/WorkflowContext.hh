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

#pragma once

#include <iostream>
#include <memory>
#include <string>

#include "common/Parameters.hh"

namespace ehunter
{

struct ContextParameters
{
    ContextParameters(HeuristicParameters heuristics)
        : heuristics(std::move(heuristics))
    {
    }
    HeuristicParameters heuristics;
};

class WorkflowContext
{
public:
    explicit WorkflowContext(std::unique_ptr<ContextParameters> paramsPtr)
    {
        if (paramsPtr == nullptr)
        {
            throw std::logic_error("Attempting to create uninitialized workflow context");
        }

        if (paramsPtr_ != nullptr)
        {
            throw std::logic_error("Attempting to redefine workflow context");
        }

        paramsPtr_ = std::move(paramsPtr);
    }

    WorkflowContext()
    {
        if (paramsPtr_ == nullptr)
        {
            throw std::runtime_error("Attempting to access uninitialized workflow context");
        }
    }

    const HeuristicParameters& heuristics() const { return paramsPtr_->heuristics; }

private:
    static std::unique_ptr<ContextParameters> paramsPtr_;
};

void initializeWorkflowContext(HeuristicParameters heuristics);

std::ostream& operator<<(std::ostream& out, const WorkflowContext& context);

}
