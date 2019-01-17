# Usage

Expansion Hunter requires an indexed BAM or a CRAM file containing aligned reads
from a PCR-free WGS sample, a FASTA file with a reference genome assembly (which
must be the same as the one used to align the reads), and a [variant catalog
file](04_VariantCatalogs.md).

Expansion Hunter outputs a VCF file and a JSON file with repeat genotypes along
with other useful information and a log file with alignments of spanning and
flanking reads and sequences of in-repeat reads. The VCF and JSON files are
largely equivalent, but the JSON file may be easier to parse programmatically.
Here is a template with the names of the required parameters.

```bash
ExpansionHunter --reads <BAM/CRAM file with aligned reads> \
                --reference <FASTA file with reference genome> \
                --variant-catalog <JSON file specifying variants to genotype> \
                --output-prefix <Prefix for the output files>
```

## Optional arguments

In addition to the required program options listed above, there are a number of
optional arguments.

* `--sex <arg>` Specifies sex of the sample; can be either `male` or `female`
  (default). This parameter only affects repeats on sex chromosomes.

* `--genome-coverage <float>` Specifies read depth on diploid chromosomes. Specifying
  read depth is required for BAM files containing a subset of alignments or CRAMs.

* `--region-extension-length <int>` Specifies how far from on/off-target regions
   to search for informative reads. Set to 1000 by default.

Note that the full list of program options with brief explanations can be
obtained by running `ExpansionHunter --help`.