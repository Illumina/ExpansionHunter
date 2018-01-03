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

#include "graphs/graph_builders.h"

#include <cmath>

using std::string;

GraphUniquePtr MakeDeletionGraph(const string& left_flank,
                                 const string& deletion,
                                 const string& right_flank) {
  GraphUniquePtr graph_ptr(new Graph(3));

  graph_ptr->SetNodeSeq(0, left_flank);
  graph_ptr->SetNodeSeq(1, deletion);
  graph_ptr->SetNodeSeq(2, right_flank);
  graph_ptr->AddEdge(0, 1);
  graph_ptr->AddEdge(0, 2);
  graph_ptr->AddEdge(1, 2);

  return graph_ptr;
}

GraphUniquePtr MakeSwapGraph(const string& left_flank, const string& deletion,
                             const string& insertion,
                             const string& right_flank) {
  GraphUniquePtr graph_ptr(new Graph(4));

  graph_ptr->SetNodeSeq(0, left_flank);
  graph_ptr->SetNodeSeq(1, deletion);
  graph_ptr->SetNodeSeq(2, insertion);
  graph_ptr->SetNodeSeq(3, right_flank);
  graph_ptr->AddEdge(0, 1);
  graph_ptr->AddEdge(0, 2);
  graph_ptr->AddEdge(1, 3);
  graph_ptr->AddEdge(2, 3);

  return graph_ptr;
}

GraphUniquePtr MakeDoubleSwapGraph(
    const string& left_flank, const string& deletion1, const string& insertion1,
    const string& middle, const string& deletion2, const string& insertion2,
    const string& right_flank) {
  GraphUniquePtr graph_ptr(new Graph(7));

  graph_ptr->SetNodeSeq(0, left_flank);
  graph_ptr->SetNodeSeq(1, deletion1);
  graph_ptr->SetNodeSeq(2, insertion1);
  graph_ptr->SetNodeSeq(3, middle);
  graph_ptr->SetNodeSeq(4, deletion2);
  graph_ptr->SetNodeSeq(5, insertion2);
  graph_ptr->SetNodeSeq(6, right_flank);
  graph_ptr->AddEdge(0, 1);
  graph_ptr->AddEdge(0, 2);
  graph_ptr->AddEdge(1, 3);
  graph_ptr->AddEdge(2, 3);
  graph_ptr->AddEdge(3, 4);
  graph_ptr->AddEdge(3, 5);
  graph_ptr->AddEdge(4, 6);
  graph_ptr->AddEdge(5, 6);

  return graph_ptr;
}

GraphUniquePtr MakeLooplessStrGraph(int32_t read_len,
                                    const std::string& left_flank,
                                    const std::string& repeat_unit,
                                    const std::string& right_flank) {
  const int32_t num_repeat_unit_nodes =
      (int)std::ceil(read_len / (double)repeat_unit.length());
  const int32_t num_nodes = num_repeat_unit_nodes + 2;  // Account for flanks

  GraphUniquePtr graph_ptr(new Graph(num_nodes));

  int32_t right_flank_node_id = num_nodes - 1;

  graph_ptr->SetNodeSeq(0, left_flank);
  graph_ptr->SetNodeSeq(right_flank_node_id, right_flank);
  graph_ptr->AddEdge(0, right_flank_node_id);

  for (int32_t node_id = 0; node_id != num_repeat_unit_nodes; ++node_id) {
    graph_ptr->SetNodeSeq(node_id + 1, repeat_unit);
    graph_ptr->AddEdge(node_id, node_id + 1);
    graph_ptr->AddEdge(node_id + 1, right_flank_node_id);
  }

  return graph_ptr;
}

GraphUniquePtr MakeStrGraph(const std::string& left_flank,
                            const std::string& repeat_unit,
                            const std::string& right_flank) {
  GraphUniquePtr graph_ptr(new Graph(3));
  graph_ptr->SetNodeSeq(0, left_flank);
  graph_ptr->SetNodeSeq(1, repeat_unit);
  graph_ptr->SetNodeSeq(2, right_flank);

  graph_ptr->AddEdge(0, 1);
  graph_ptr->AddEdge(0, 2);
  graph_ptr->AddEdge(1, 1);
  graph_ptr->AddEdge(1, 2);

  return graph_ptr;
}
