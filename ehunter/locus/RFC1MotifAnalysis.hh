//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Chris Saunders <csaunders@illumina.com>
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

#include "locus/AlignmentBuffer.hh"
#include "locus/LocusFindings.hh"

namespace ehunter
{

/// \brief Analyze the RFC1 locus motif pattern with respect to motif expansions associated with CANVAS
///
/// RFC1 locus call information are added to LocusFindings for inclusion in EH json output.
///
/// \param[in] alignmentBuffer
///
/// \param[in,out] locusFindings Findings from conventional repeat expansion analysis of RFC1. Additional RFC1 motif
/// analysis is added to this object for reporting downstream in the EH json output.
///
void runRFC1MotifAnalysis(const locus::AlignmentBuffer& alignmentBuffer, LocusFindings& locusFindings);

}
