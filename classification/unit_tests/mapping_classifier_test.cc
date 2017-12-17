//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
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

#include "classification/mapping_classifier.h"

#include <memory>
#include <string>

#include "gtest/gtest.h"

#include "graphs/graph.h"
#include "graphs/graph_builders.h"
#include "graphs/graph_mapping.h"
#include "graphs/graph_mapping_operations.h"

using std::string;

TEST(MappingClassificaton, SpanningRead_Classified) {
  Graph graph = MakeStrGraph("AAAACC", "CCG", "ATTT");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);
  //                            FFRRRRRRFF
  const string spanning_read = "CCCCGCCGAT";
  GraphMapping spanning_mapping =
      DecodeFromString(4, "0[2M]1[3M]1[3M]2[2M]", spanning_read, graph);

  StrMappingClassifier mapping_classifier(0, 1, 2);
  ASSERT_EQ(ReadClass::kSpansRepeat,
            mapping_classifier.Classify(spanning_mapping));
}