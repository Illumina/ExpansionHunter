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
#include "genotyping/AlleleChecker.hh"
#include "genotyping/RepeatGenotype.hh"
#include "genotyping/SmallVariantGenotype.hh"

namespace ehunter
{

class StrFindings;
class SmallVariantFindings;
class CnvVariantFindings;

struct VariantFindingsVisitor
{
    virtual void visit(StrFindings& strFindings) = 0;
    virtual void visit(SmallVariantFindings& smallVariantFindings) = 0;
    virtual void visit(CnvVariantFindings& cnvVariantFindings) = 0;
};

class VariantFindings
{
public:
    VariantFindings(std::string variantId)
        : variantId_(std::move(variantId))
    {
    }

    virtual ~VariantFindings() = default;
    virtual void accept(VariantFindingsVisitor& visitor) = 0;
    const std::string& variantId() const { return variantId_; }

protected:
    std::string variantId_;
};

class StrFindings : public VariantFindings
{
public:
    StrFindings(
        std::string variantId, CountTable countsOfSpanningReads, CountTable countsOfFlankingReads,
        CountTable countsOfInrepeatReads, boost::optional<RepeatGenotype> optionalGenotype)
        : VariantFindings(std::move(variantId))
        , countsOfSpanningReads_(std::move(countsOfSpanningReads))
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

class SmallVariantFindings : public VariantFindings
{
public:
    SmallVariantFindings(
        std::string variantId, int numRefReads, int numAltReads, AlleleCheckSummary refAlleleStatus,
        AlleleCheckSummary altAlleleStatus, boost::optional<SmallVariantGenotype> optionalGenotype)
        : VariantFindings(std::move(variantId))
        , numRefReads_(numRefReads)
        , numAltReads_(numAltReads)
        , refAlleleStatus_(refAlleleStatus)
        , altAlleleStatus_(altAlleleStatus)
        , optionalGenotype_(std::move(optionalGenotype))
    {
    }

    ~SmallVariantFindings() override = default;
    void accept(VariantFindingsVisitor& visitorPtr) override { visitorPtr.visit(*this); }

    int numRefReads() const { return numRefReads_; }
    int numAltReads() const { return numAltReads_; }
    const boost::optional<SmallVariantGenotype>& optionalGenotype() const { return optionalGenotype_; }

    AlleleCheckSummary refAllelePresenceStatus() const { return refAlleleStatus_; }
    AlleleCheckSummary altAllelePresenceStatus() const { return altAlleleStatus_; }

private:
    int numRefReads_;
    int numAltReads_;
    AlleleCheckSummary refAlleleStatus_;
    AlleleCheckSummary altAlleleStatus_;
    boost::optional<SmallVariantGenotype> optionalGenotype_;
};

class CnvVariantFindings : public VariantFindings
{
public:
    CnvVariantFindings(std::string variantId, boost::optional<int> copyNumberCall)
        : VariantFindings(std::move(variantId))
        , copyNumberCall_(copyNumberCall)
    {
    }

    ~CnvVariantFindings() override = default;
    void accept(VariantFindingsVisitor& visitorPtr) override { visitorPtr.visit(*this); }

    boost::optional<int> copyNumberCall() const { return copyNumberCall_; }

private:
    boost::optional<int> copyNumberCall_;
};
}
