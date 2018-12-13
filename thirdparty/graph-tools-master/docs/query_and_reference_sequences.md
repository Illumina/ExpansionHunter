# Query and reference sequences

## Introduction

A large part of the GraphTools library is concerned with alignment of query sequences (typically short reads) to a
reference represented by a sequence graph. This section explains the distinction between query sequences and reference
sequences that make up graphs.

A **query sequence** must consist of characters representing core nucleotides A, C, G, T and a missing nucleotide symbol
N. The uppercase and lowercase core nucleotide characters correspond to high and low quality base calls respectively.

Here is an example of a query sequence with many low quality base calls.

```C++
CCGACCACGCCCCGGCCCCcGCCCCGGCCCCcaGCGcgCGaCcCCtGaGgTcccGgGctTGCcaCaGgCcggcGgtGttTCCCcCCttgttttTTtCtg
```

A **reference sequence** -- a sequence of a node in a graph -- must consist of characters representing cores nucleotides
(A, C, G, T) or degenerate nucleotides (that correspond to multiple nucleotides). The following table lists the
degenerate nucleotides as defined by the International Union of Pure and Applied Chemistry.

| Degenerate nucleotide | Corresponding core nucleotides |
|-----------------------|--------------------------------|
| R                     | A, G                           |
| Y                     | C, T                           |
| K                     | G, T                           |
| M                     | A, C                           |
| S                     | C, G                           |
| W                     | A, T                           |
| B                     | C, G, T                        |
| D                     | A, G, T                        |
| H                     | A, C, T                        |
| V                     | A, C, G                        |
| N                     | A, C, G, T                     |

Here is an example of a reference sequence representing a trinucleotide repeat consisting of 12 Alanine codons.

```C++
GCNGCNGCNGCNGCNGCNGCNGCNGCNGCNGCNGCN
```

## Matching query and reference sequences

The rules for matching the characters of query and reference sequences are given in the table above. Additionally, we
require that the missing nucleotide character N in the query sequences does not match any other characters. (Note the
distinction with reference N that matches A, C, G, and T.)

To accommodate the needs of the existing clients of the GraphTools library, we allow query Ns within alignments to be a
part of mismatch operations or "missing nucleotide" operations. See the section on
[representation of alignments](representation_of_alignments.md) for more details.

The basic machinery for matching query and reference sequences is provided in `BaseMatching.*` and
`SequenceOperations.*` files inside `graphutils`.
