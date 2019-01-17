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

#include <memory>
#include <unordered_map>

#include <boost/optional.hpp>

#include "common/CountTable.hh"
#include "genotyping/AllelePresenceChecker.hh"
#include "genotyping/RepeatGenotype.hh"
#include "genotyping/SmallVariantGenotype.hh"

namespace ehunter
{

class VariantFindings;
class RepeatFindings;
class SmallVariantFindings;

struct VariantFindingsVisitor
{
    virtual void visit(const RepeatFindings* findingsPtr) = 0;
    virtual void visit(const SmallVariantFindings* findingsPtr) = 0;
};

class VariantFindings
{
public:
    virtual ~VariantFindings() = default;
    virtual void accept(VariantFindingsVisitor* visitorPtr) = 0;
};

using RegionFindings = std::unordered_map<std::string, std::unique_ptr<VariantFindings>>;
using SampleFindings = std::unordered_map<std::string, RegionFindings>;

class RepeatFindings : public VariantFindings
{
public:
    RepeatFindings(
        CountTable countsOfSpanningReads, CountTable countsOfFlankingReads, CountTable countsOfInrepeatReads,
        boost::optional<RepeatGenotype> optionalGenotype)
        : countsOfSpanningReads_(std::move(countsOfSpanningReads))
        , countsOfFlankingReads_(std::move(countsOfFlankingReads))
        , countsOfInrepeatReads_(std::move(countsOfInrepeatReads))
        , optionalGenotype_(std::move(optionalGenotype))
    {
    }

    ~RepeatFindings() override = default;
    void accept(VariantFindingsVisitor* visitorPtr) override { visitorPtr->visit(this); }

    const CountTable& countsOfSpanningReads() const { return countsOfSpanningReads_; }
    const CountTable& countsOfFlankingReads() const { return countsOfFlankingReads_; }
    const CountTable& countsOfInrepeatReads() const { return countsOfInrepeatReads_; }
    const boost::optional<RepeatGenotype>& optionalGenotype() const { return optionalGenotype_; }

    bool operator==(const RepeatFindings& other) const
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
        int numRefReads, int numAltReads,
        AllelePresenceStatus refAlleleStatus, AllelePresenceStatus altAlleleStatus,
        boost::optional<SmallVariantGenotype> optionalGenotype)
        : numRefReads_(numRefReads)
        , numAltReads_(numAltReads)
        , refAlleleStatus_(refAlleleStatus)
        , altAlleleStatus_(altAlleleStatus)
        , optionalGenotype_(std::move(optionalGenotype))
    {
    }

    ~SmallVariantFindings() override = default;
    void accept(VariantFindingsVisitor* visitorPtr) override { visitorPtr->visit(this); }

    int numRefReads() const { return numRefReads_; }
    int numAltReads() const { return numAltReads_; }
    const boost::optional<SmallVariantGenotype>& optionalGenotype() const { return optionalGenotype_; }

    AllelePresenceStatus refAllelePresenceStatus() const { return refAlleleStatus_; }
    AllelePresenceStatus altAllelePresenceStatus() const { return altAlleleStatus_; }

private:
    int numRefReads_;
    int numAltReads_;
    AllelePresenceStatus refAlleleStatus_;
    AllelePresenceStatus altAlleleStatus_;
    boost::optional<SmallVariantGenotype> optionalGenotype_;
};

std::ostream& operator<<(std::ostream& out, const RepeatFindings& repeatFindings);

}
