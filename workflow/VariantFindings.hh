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
#include "common/GenomicRegion.hh"
#include "genotyping/AlleleChecker.hh"
#include "genotyping/RepeatGenotype.hh"
#include "genotyping/SmallVariantGenotype.hh"

namespace ehunter
{

class StrFindings;
class SmallVariantFindings;
class CnvVariantFindings;
class ParalogSmallVariantFindings;

struct VariantFindingsVisitor
{
    virtual void visit(const StrFindings& findings) = 0;
    virtual void visit(const SmallVariantFindings& findings) = 0;
    virtual void visit(const CnvVariantFindings& findings) = 0;
    virtual void visit(const ParalogSmallVariantFindings& findings) = 0;
};

class VariantFindings
{
public:
    explicit VariantFindings(std::string variantId)
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

class ParalogSmallVariantFindings : public VariantFindings
{
public:
    ParalogSmallVariantFindings(
        std::string variantId, int numGeneAReads, int numGeneBReads, boost::optional<std::pair<int, double>> copyNumber)
        : VariantFindings(std::move(variantId))
        , numGeneAReads_(numGeneAReads)
        , numGeneBReads_(numGeneBReads)
        , copyNumber_(std::move(copyNumber))
    {
    }

    ~ParalogSmallVariantFindings() override = default;
    void accept(VariantFindingsVisitor& visitorPtr) override { visitorPtr.visit(*this); }

    int numGeneAReads() const { return numGeneAReads_; }
    int numGeneBReads() const { return numGeneBReads_; }
    const boost::optional<std::pair<int, double>>& copyNumber() const { return copyNumber_; }

private:
    int numGeneAReads_;
    int numGeneBReads_;
    boost::optional<std::pair<int, double>> copyNumber_;
};

class CnvVariantFindings : public VariantFindings
{
public:
    CnvVariantFindings(std::string variantId, boost::optional<int> absoluteCopyNumber, boost::optional<int> copyNumberChange)
        : VariantFindings(std::move(variantId))
        , absoluteCopyNumber_(absoluteCopyNumber)
        , copyNumberChange_(copyNumberChange)
    {
    }

    ~CnvVariantFindings() override = default;
    void accept(VariantFindingsVisitor& visitorPtr) override { visitorPtr.visit(*this); }

    boost::optional<int> absoluteCopyNumber() const { return absoluteCopyNumber_; }
    boost::optional<int> copyNumberChange() const { return copyNumberChange_; }

private:
    boost::optional<int> absoluteCopyNumber_;
    boost::optional<int> copyNumberChange_;
};
}
