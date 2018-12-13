# Graph alignment

## Alignment methods

## Alignment representation

## BAM Output

Graph-tools supports outputting alignments in BAM format (through `graphIO/BamWriter`).
Currently this does not produce full BAM alignments, but rather 'unaligned' BAM records 
that are placed at the right start position in the (linear) reference genome and include
the actual alignment against the alignment encoded as a string in a custom tag.

For this to work a projection of the graph onto a linear reference genome is needed (class `GraphReferenceMapping`).
This assigns a reference coordinate interval to graph nodes and hence allows projecting positions (and paths) in the graph to
positions in the linear reference that the BAM is based on.

### BAM encoding for graph alignments

The BAM record for read `R` with graph-alignment `g`:
  * Qname, Seq, Qual describe the read (matching fastq input)
  * Chromosome and Position are set to the first projected reference position of the graph-alignment path (`g.path`)
    1. Find first node `n` in `g.path` that has a reference mapping
    2. Set the BAM read position to the projected position of the first base of `n` + the offset of `g.path` in `n`.
  * Mate Chrom and Position undefined
  * Read is `unmapped`: No mapQ, CIGAR
  * Flags
    * read paired, first/second in pair, failedQC: Set according to fastq/input BAM
    * read unmapped and mate unmapped is set
    * read reverse strand is set depending on mapping to graph in forward or reverse direction
    * all other flags are unset
  * The actual graph alignment is encoded in custom tags
    * XG: string - Graph CIGAR
      * Combination of Alignment path and CIGAR against each node.
      * E.g. 0[Ref start: 2, 2M]1[Ref start: 0, 1M]3[Ref start: 0, 3M]
        * Alignment starts at position 2 on node 0 with 2 matches
        * Continues at Node1 with 1 match and node3 with 3 matches
