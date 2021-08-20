# Usage

Expansion Hunter requires an indexed BAM or a CRAM file containing aligned reads
from a PCR-free WGS sample, a FASTA file with a reference genome assembly (which
must be the same as the one used to align the reads), and a [variant catalog
file](04_VariantCatalogFiles.md).

Expansion Hunter outputs a VCF file and a JSON file with variant genotypes and
other useful information along with a BAMlet containing alignments of reads that
overlap or located in close proximity to each variant. The VCF and JSON files
are largely equivalent, but the JSON file may be easier to parse
programmatically. Here is a template with the names of the required parameters.

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

* `--region-extension-length <int>` Specifies how far from on/off-target regions
   to search for informative reads. Set to 1000 by default.

* `--analysis-mode <mode>` Specify analysis mode, which can be either `seeking` or
  `streaming`. The default mode is `seeking`. See further description of analysis
   modes below.

* `--threads <int>` Specifies how many threads to can be used accelerate analysis
   of large variant catalogs. Set to 1 by default. Typically seeking mode can
   benefit from relatively high thread counts, while for streaming mode
   there is limited benefit beyond about 16 threads.

Note that the full list of program options with brief explanations can be
obtained by running `ExpansionHunter --help`.

### Analysis modes

#### Seeking mode

In seeking mode, alignment file indexing is used to seek specific read sets for the
analysis of each variant. Seeking mode is recommended for analysis of small catalogs.
This mode requires that the input BAM or CRAM file is already sorted and indexed.

#### Streaming mode

In streaming mode, the alignment file is read in a single pass and all variants are
analyzed during this reading operation. Streaming mode is recommended for the analysis
of large catalogs, but does require more memory as a funciton of catalog size. This mode
does not require that the BAM or CRAM file is sorted or indexed.
