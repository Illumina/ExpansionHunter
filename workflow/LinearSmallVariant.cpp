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

#include "workflow/LinearSmallVariant.hh"

#include "workflow/LinearModel.hh"

using std::shared_ptr;
using std::vector;

namespace ehunter
{

LinearSmallVariant::LinearSmallVariant(
    std::shared_ptr<LinearModel> model, SmallVariantLocations locations, SmallVariantBases bases, boost::optional<int> mapqCutoff)
    : model_(std::move(model))
    , locations_(std::move(locations))
    , bases_(std::move(bases))
    , mapqCutoff_(std::move(mapqCutoff))
{
}

shared_ptr<RegionModel> LinearSmallVariant::model() { return model_; }

static boost::optional<Base> decodeBase(const char base)
{
    if (base == 'A')
    {
        return Base::kA;
    }
    else if (base == 'C')
    {
        return Base::kC;
    }
    else if (base == 'T')
    {
        return Base::kT;
    }
    else if (base == 'G')
    {
        return Base::kG;
    }

    return boost::optional<Base>();
}

static boost::optional<Base> getBaseOnRead(const MappedRead& read, int position)
{
    std::vector<std::pair<char, int>> cigarOps = read.cigarOp();
    
    int positionOnReference = position;
    int positionOnQuery = 0;
    // int referenceOnly = 0;

    for (auto cigarOp : cigarOps)
    {
        if (cigarOp.first == 'S' || cigarOp.first == 'I')
        {
            positionOnQuery += cigarOp.second;
        }
        if (cigarOp.first == 'M' || cigarOp.first == '=' || cigarOp.first == 'X')
        {
            positionOnQuery += cigarOp.second;
            positionOnReference += cigarOp.second;
        }
        if (cigarOp.first == 'D')
        {
            positionOnReference += cigarOp.second;
        }

        if (positionOnReference <= position) continue;

        if (cigarOp.first == 'M' || cigarOp.first == '=' || cigarOp.first == 'X')
        {
            int NumberOfBasesPastPosition = positionOnReference - position;
            int basePositionOnQuery = positionOnQuery - NumberOfBasesPastPosition;
            if (basePositionOnQuery >= (int)read.sequence().size() || basePositionOnQuery < 0)
            {
                throw std::logic_error("Position past read end.");
            }
            char readBase = read.sequence()[basePositionOnQuery];
            std::cout << read.readId() << " " << basePositionOnQuery << " " << position << " " << readBase << "\n"; 
            return decodeBase(readBase);
        }
        return boost::optional<Base>();
    }
    return boost::optional<Base>();
}

void LinearSmallVariant::summarize(const MappedRead& read)
{
    if (mapqCutoff_)
    {
        if (read.mapq() >= *mapqCutoff_)
        {
            int posA = locations_.geneALocation.start();
            int posB = locations_.geneBLocation.start();
            
            boost::optional<Base> variantBase;
            if (read.pos() < posA && read.approximateEnd() + 100 > posA)
            {
                std::cout << read.readId() << " " << read.pos() << " " << posA << " 1 \n";
                variantBase = getBaseOnRead(read, posA);
            }
            else if (read.pos() < posB && read.approximateEnd() + 100 > posB)
            {
                std::cout << read.readId() << " " << read.pos() << " " << posB << " 2 \n";
                variantBase = getBaseOnRead(read, posB);
            }
            if (variantBase)
            {
                if (*variantBase == bases_.geneABase)
                {
                    ++numGeneAReads_;
                }
                else if (*variantBase == bases_.geneBBase)
                {
                    ++numGeneBReads_;
                }
            }
        }
    }
}

void LinearSmallVariant::summarize(const MappedRead& read, const MappedRead& mate)
{
    summarize(read);
    summarize(mate);
}
}