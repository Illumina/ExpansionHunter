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

#include <list>
#include <memory>
#include <string>

#include "gtest/gtest.h"

#include "graphs/graph.h"
#include "graphs/graph_builders.h"
#include "graphs/graph_mapping.h"
#include "graphs/graph_mapping_operations.h"

using std::list;
using std::string;

TEST(MappingClassificaton, SpanningMapping_Classified) {
  Graph graph = MakeStrGraph("AAAACC", "CCG", "ATTT");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);
  StrMappingClassifier mapping_classifier(0, 1, 2);

  {  //                  FFRRRRRRFF
    const string read = "CCCCGCCGAT";
    GraphMapping mapping =
        DecodeFromString(4, "0[2M]1[3M]1[3M]2[2M]", read, graph);

    EXPECT_EQ(MappingType::kSpansRepeat, mapping_classifier.Classify(mapping));
  }

  {  //                  FFFF
    const string read = "CCAT";
    GraphMapping mapping = DecodeFromString(4, "0[2M]2[2M]", read, graph);

    StrMappingClassifier mapping_classifier(0, 1, 2);
    EXPECT_EQ(MappingType::kSpansRepeat, mapping_classifier.Classify(mapping));
  }
}

TEST(MappingClassificaton, FlankingMapping_Classified) {
  Graph graph = MakeStrGraph("AAAACC", "CCG", "ATTT");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);
  StrMappingClassifier mapping_classifier(0, 1, 2);

  {  //                  FFFFRRR
    const string read = "AACCCCG";
    GraphMapping mapping = DecodeFromString(2, "0[4M]1[3M]", read, graph);

    EXPECT_EQ(MappingType::kFlanksRepeat, mapping_classifier.Classify(mapping));
  }

  {  //                  RRRFFF
    const string read = "CCGATT";
    GraphMapping mapping = DecodeFromString(0, "1[3M]2[3M]", read, graph);

    StrMappingClassifier mapping_classifier(0, 1, 2);
    EXPECT_EQ(MappingType::kFlanksRepeat, mapping_classifier.Classify(mapping));
  }
}

TEST(MappingClassificaton, RepeatMapping_Classified) {
  Graph graph = MakeStrGraph("AAAACC", "CCG", "ATTT");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);
  StrMappingClassifier mapping_classifier(0, 1, 2);

  {  //                  RRRRRRRR
    const string read = "CCGCCGCC";
    GraphMapping mapping = DecodeFromString(0, "1[3M]1[3M]1[2M]", read, graph);

    EXPECT_EQ(MappingType::kInsideRepeat, mapping_classifier.Classify(mapping));
  }

  {  //                  RRRRRRRR
    const string read = "CGCCGCCG";
    GraphMapping mapping = DecodeFromString(1, "1[2M]1[3M]1[3M]", read, graph);

    StrMappingClassifier mapping_classifier(0, 1, 2);
    EXPECT_EQ(MappingType::kInsideRepeat, mapping_classifier.Classify(mapping));
  }
}

TEST(MappingClassificaton, OutsideRepeatMapping_Classified) {
  Graph graph = MakeStrGraph("AAAACC", "CCG", "ATTT");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);
  StrMappingClassifier mapping_classifier(0, 1, 2);

  {  //                  FFFFF
    const string read = "AAAAC";
    GraphMapping mapping = DecodeFromString(0, "0[5M]", read, graph);

    EXPECT_EQ(MappingType::kOutsideRepeat,
              mapping_classifier.Classify(mapping));
  }

  {  //                  FFF
    const string read = "TTT";
    GraphMapping mapping = DecodeFromString(1, "2[3M]", read, graph);

    StrMappingClassifier mapping_classifier(0, 1, 2);
    EXPECT_EQ(MappingType::kOutsideRepeat,
              mapping_classifier.Classify(mapping));
  }
}

TEST(ObtainingCanonicalMapping, FlankingAndInrepeatRead_ClassifiedAsInrepeat) {
  Graph graph = MakeStrGraph("AAAACG", "CCG", "ATTT");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);
  StrMappingClassifier mapping_classifier(0, 1, 2);

  //                   FFFFFFFF
  const string read = "CGCCGCCG";
  const GraphMapping flanking_mapping =
      DecodeFromString(4, "0[2M]1[3M]1[3M]", read, graph);
  const GraphMapping irr_mapping =
      DecodeFromString(1, "1[2M]1[3M]1[3M]", read, graph);

  const list<GraphMapping> mappings = {flanking_mapping, irr_mapping};

  EXPECT_EQ(irr_mapping, mapping_classifier.GetCanonicalMapping(mappings));
}

TEST(ObtainingCanonicalMapping, FlankingAndSpanningRead_ClassifiedAsFlanking) {
  Graph graph = MakeStrGraph("AAAACG", "CCG", "ATTT");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);
  StrMappingClassifier mapping_classifier(0, 1, 2);

  //                   FFFFFFFFFF
  const string read = "CGCCGCCGAT";
  const GraphMapping spanning_mapping =
      DecodeFromString(4, "0[2M]1[3M]1[3M]2[2M]", read, graph);
  const GraphMapping flanking_mapping =
      DecodeFromString(1, "1[2M]1[3M]1[3M]2[2M]", read, graph);

  const list<GraphMapping> mappings = {spanning_mapping, flanking_mapping};

  EXPECT_EQ(flanking_mapping, mapping_classifier.GetCanonicalMapping(mappings));
}