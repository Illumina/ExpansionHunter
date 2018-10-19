# Graph model

## Definitions
The sequence graph is made up of nodes and directed edges. Each node holds a (non-empty) sequence. Edges connect
nodes and do not hold sequences themselves. A path is a sequence of nodes, where each node is connected to its successor with a directed edge. A path defines a sequence in the graph.

## Graph class
### Files
- `include/graphcore/Graph.hh`
- `tests`
    - `GraphOperationsTest.cpp`
    - `GraphCoordinatesTest.cpp`

The graph class holds nodes, edges and edge labels (used for path families).
Nodes are simple structs with a name and sequence. Within a graph each node has a unique NodeID, numbered from 0 to
graph.numNodes.
Edges are defined as pairs of NodeIDs and are stored directly within the graph using a node adjacency list.
Every edge can optionally be labeled with one or more labels.

Graphs are allowed to contain self-loops (edges from a node to itself). Excluding self-loop edges the graph has to be a DAG. A path
may pass through such an loop-edge multiple times, defining a repeated sequence.

## Paths
See [path](paths.md) and [path family](path_families.md) documentation.

## Reference Mapping
A graph reference mapping (`include/Graphcore/GraphReferenceMapping.hh`) is used to define a projection of the graph onto linear reference genome coordinates. A subset of nodes
are assigned one (or possibly multiple) reference intervals. For reference nodes this would be the interval their sequence comes from
for alt. allele nodes the 'replaced' reference sequence (or position of insertion).

These reference mappings can be stored in JSON, either together with the graph or separately. They are used to generate BAM output of alignments.

## JSON model
Graphs can be serialized in JSON. The format is
- "nodes": array of node objects containing
    - "name": string (required; not-empty)  - Unique name for this node
    - "sequence": string; not-empty - Nucleotide sequence of this node
        - [Extended alphabet](query_and_reference_sequences.md)
        - Either "sequence" or "reference" is required
    - "reference": string - Region in the reference genome
        - Sets the sequence of the node to match the reference region if no explicit "sequence" is given
        - Region is given as a 0-based half-open interval (<chrom>:<Start>-<End (excluded)>)
        - Also defines this reference region as the mapping position of the node for GraphReferenceMapping
- "edges": array of edge objects containing
    - "from": string - name of the start node
    - "to": string - name of the end node
    - "labels": array of string - edge labels applied to this edge
- "paths": array of path objects (optional) containing
    - path_id: string - name of the path
    - nodes: array of string - name of nodes in the path (in order)
- "graph_id": string (optional) - name for the graph
- "reference_genome": string - path to a fasta file with the reference genome sequence
    - required if using "reference" sequences in node definitions
