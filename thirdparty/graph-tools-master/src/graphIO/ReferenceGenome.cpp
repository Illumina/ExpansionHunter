// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "graphIO/ReferenceGenome.hh"

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using std::string;

namespace graphIO
{

RefGenome::RefGenome(string const& fastaPath)
    : fastaPath_(fastaPath)
    , fai_(fai_load(fastaPath.c_str()), fai_destroy)
{
}

string RefGenome::extractSeq(ReferenceInterval const& interval) const
{
    int len;
    // pass end - 1 since htslib includes the last base, but our interval object excludes it
    std::unique_ptr<char[]> refTmp(
        faidx_fetch_seq(fai_.get(), interval.contig.c_str(), interval.start, interval.end - 1, &len));
    if (!refTmp || len == -1 || len == -2 || len != interval.length())
    {
        throw std::runtime_error((boost::format("ERROR: can't extract %1% from %2%") % interval % fastaPath_).str());
    }
    string seq(refTmp.get(), len);
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    return seq;
}
}