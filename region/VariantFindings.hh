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

#include <boost/optional.hpp>

#include "common/CountTable.hh"
#include "genotyping/RepeatGenotype.hh"

namespace ehunter
{

class StrFindings;

struct VariantFindingsVisitor
{
    virtual void visit(StrFindings& strFindings) = 0;
};

struct VariantFindings
{
    virtual ~VariantFindings() = default;
    virtual void accept(VariantFindingsVisitor& visitor) = 0;
};

class StrFindings : public VariantFindings
{
public:
    StrFindings(
        CountTable countsOfSpanningReads, CountTable countsOfFlankingReads, CountTable countsOfInrepeatReads,
        boost::optional<RepeatGenotype> optionalGenotype)
        : countsOfSpanningReads_(std::move(countsOfSpanningReads))
        , countsOfFlankingReads_(std::move(countsOfFlankingReads))
        , countsOfInrepeatReads_(std::move(countsOfInrepeatReads))
        , optionalGenotype_(std::move(optionalGenotype))
    {
    }

    ~StrFindings() override = default;
    void accept(VariantFindingsVisitor& visitor) override { visitor.visit(*this); }

    const CountTable& countsOfSpanningReads() const { return countsOfSpanningReads_; }
    const CountTable& countsOfFlankingReads() const { return countsOfFlankingReads_; }
    const CountTable& countsOfInrepeatReads() const { return countsOfInrepeatReads_; }
    const boost::optional<RepeatGenotype>& optionalGenotype() const { return optionalGenotype_; }

    bool operator==(const StrFindings& other) const
    {
        return countsOfSpanningReads_ == other.countsOfSpanningReads_
            && countsOfFlankingReads_ == other.countsOfFlankingReads_
            && countsOfInrepeatReads_ == other.countsOfInrepeatReads_ && optionalGenotype_ == other.optionalGenotype_;
    }

private:
    CountTable countsOfSpanningReads_;
    CountTable countsOfFlankingReads_;
    CountTable countsOfInrepeatReads_;
    boost::optional<RepeatGenotype> optionalGenotype_;
};

}
