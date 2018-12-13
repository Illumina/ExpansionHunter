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

#include "common/Parameters.hh"

#include <iomanip>
#include <ios>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

using std::string;
using std::vector;

namespace ehunter {


namespace po = boost::program_options;
namespace fs = boost::filesystem;

Outputs::Outputs(const string& vcfPath, const string& jsonPath, const string& logPath)
{
    vcf_.open(vcfPath.c_str());
    if (!vcf_.is_open())
    {
        throw std::runtime_error("Failed to open " + vcfPath + " for writing: " + strerror(errno));
    }

    json_.open(jsonPath.c_str());
    if (!json_.is_open())
    {
        throw std::runtime_error("Failed to open " + jsonPath + " for writing: " + strerror(errno));
    }

    log_.open(logPath.c_str());
    if (!log_.is_open())
    {
        throw std::runtime_error("Failed to open " + logPath + " for writing: " + strerror(errno));
    }
}


}
