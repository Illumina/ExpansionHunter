# Usage

Expansion Hunter requires an indexed BAM or a CRAM file containing aligned reads from a PCR-free WGS sample, a FASTA
file with a reference genome assembly (which must be identical to the one used for aligning the reads), and repeat
specification files ([see here](Inputs)).

Expansion Hunter outputs a VCF file and a JSON file with repeat genotypes along with other useful information and a log
file with alignments of spanning and flanking reads and sequences of in-repeat reads. The VCF and JSON files are largely
equivalent, but the JSON file may be easier to parse programmatically. Here is a template that provides the names of the
required parameters.

```bash
ExpansionHunter --bam <BAM/CRAM file with aligned reads> \
                --ref-fasta <FASTA file with reference genome> \
                --repeat-specs <Directory with repeat specification files> \
                --vcf <Output VCF file> \
                --json <Output JSON file> \
                --log <Output file with alignments of relevant reads>
```

Optional arguments
------------------

In addition to the required program options listed above, there are a number of optional arguments.

| Program option  | Description  | Default |
|-----------------|--------------|:-------:|
| --sex arg | Specifies sex of the sample; can be either male or female*. | female |
| --skip-unaligned | If set, requires the program to skip the search for IRRs in unaligned reads**. | -- |
| --read-depth arg | Specifies read depth so that Expansion Hunter can skip computing it***. | -- |
| --min-score arg | The minimum weighted matching score (can be any positive number less than 1). | 0.90 |
| --min-baseq arg | The minimum base quality of a high-confidence base call. | 20 |
| --min-anchor-mapq arg | The minimum MAPQ of an in-repeat read anchor | 60 |
| --region-extension-length arg | Specifies how far from on/off-target regions to search for informative reads. | 1000 |

\* `sex` parameter only affects genotyping of sex chromosomes.

\** Skipping unaligned reads can results in a considerable speed up but can also cause repeat sizes to get
underestimated.

\*** Specifying read depth is required when analyzing BAM files that contain a subset of read alignments.

Note that the full list of program options with brief explanations can be obtained by running `ExpansionHunter --help`.