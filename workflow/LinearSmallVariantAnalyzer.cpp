//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Xiao Chen <xchen2@illumina.com>,
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

#include "workflow/LinearSmallVariantAnalyzer.hh"
#include "genotyping/SmallVariantCopyNumberGenotyper.hh"

#include "common/Common.hh"

using std::shared_ptr;
using std::vector;

namespace ehunter
{
LinearSmallVariantAnalyzer::LinearSmallVariantAnalyzer(
    std::string variantId, std::shared_ptr<LinearSmallVariant> linearSmallVariant)
    : variantId_(variantId)
    , linearSmallVariant_(std::move(linearSmallVariant))
{
}

vector<shared_ptr<Feature>> LinearSmallVariantAnalyzer::features() { return { linearSmallVariant_ }; }

ParalogSmallVariantFindings LinearSmallVariantAnalyzer::analyze(boost::optional<int> totalCopyNumber) const
{
    const int numGeneAReads = linearSmallVariant_->numGeneAReads();
    const int numGeneBReads = linearSmallVariant_->numGeneBReads();
    if (totalCopyNumber)
    {
        SmallVariantCopyNumberGenotyper genotyper(*totalCopyNumber);
        auto smallVariantGenotype = genotyper.genotype(numGeneAReads, numGeneBReads, 1);
        return ParalogSmallVariantFindings(variantId_, numGeneAReads, numGeneBReads, smallVariantGenotype);
    }

    return ParalogSmallVariantFindings(variantId_, numGeneAReads, numGeneBReads, boost::none); 
    
}
}
