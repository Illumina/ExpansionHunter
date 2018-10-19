# A library for working with sequence graphs

## Introduction

Graph-tools library aims to provide data models and methods for working with sequence graphs and alignments to them.

See docs/development.md for build instructions

## Contents

### Graph definition
- A data model for graphs in memory and as JSON: [graphs.md](docs/graphs.md)
    - The graphs are DAGs, with the addition of support for self-loops on individual nodes
- Paths through graphs to define sequences, e.g. for alignments: [paths.md](docs/paths.md)
- Path Families, subsets of edges of the graph, to define haplotypes: [path_families.md](docs/path_families.md)
- Graph-to-reference mappings: Translating positions in the graph to positions in a linear reference (e.g. genome)
- Nucleotide sequence definition:[query_and_reference_sequences.md](docs/query_and_reference_sequences.md)
    - Extended alphabet for degenerate bases and base quality

### Alignment
- Methods to perform both gapped and gapless alignment of sequences to graphs,
- A rich set of heuristics for filtering alignments.

### GraphIO
The separate GraphIO library contains methods to read and write
    - Graphs from JSON
    - (Linear) Reference sequences from fasta
    - Graph-alignments from BAM
It depends on htslib.