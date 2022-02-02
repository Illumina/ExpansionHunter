# Usage

Expansion Hunter requires the following inputs:
1. A BAM or a CRAM file containing aligned reads from a PCR-free WGS sample.
    1. The BAM or CRAM file must be sorted and indexed if using the seeking [analysis mode](#analysis-modes).
    2. The BAM or CRAM file may be a local filesystem path or [URL](#url-support).
4. A FASTA file with a reference genome assembly (which must be the same as the one used to align the reads)
5. A [variant catalog file](04_VariantCatalogFiles.md).

Expansion Hunter outputs a VCF file and a JSON file with variant genotypes and
other useful information along with a BAMlet containing alignments of reads that
overlap or located in close proximity to each variant. The VCF and JSON files
are largely equivalent, but the JSON file may be easier to parse
programmatically. Here is a template with the names of the required parameters.

```bash
ExpansionHunter --reads <aligned reads BAM/CRAM file/URL> \
                --reference <reference genome FASTA file> \
                --variant-catalog <JSON file specifying variants to genotype> \
                --output-prefix <Prefix for the output files>
```

## Optional arguments

In addition to the required program options listed above, there are a number of
optional arguments.

* `--sex <arg>` Specifies sex of the sample; can be either `male` or `female`
  (default). This parameter only affects repeats on sex chromosomes.

* `--threads <int>` Specifies how many threads to can be used accelerate analysis
   of large variant catalogs. Set to 1 by default. Typically seeking mode can
   benefit from relatively high thread counts, while for streaming mode
   there is limited benefit beyond about 16 threads.

* `--min-locus-coverage <int>` Specifies minimum read coverage depth at loci
   on diploid chromosomes required to attempt genotyping. Automatically reduced
   to half for loci on haploid chromosomes. The locus will be skipped if the
   coverage falls below this value. Set to 10 by default.

* `--region-extension-length <int>` Specifies how far from on/off-target regions
   to search for informative reads. Set to 1000 by default.

* `--analysis-mode <mode>` Specify analysis mode, which can be either `seeking` or
  `streaming`. The default mode is `seeking`. See further description of analysis
   modes below.


Note that the full list of program options with brief explanations can be
obtained by running `ExpansionHunter --help`.

### URL support

The aligned reads input BAM or CRAM file may be a local filesystem path or URL.
Supported protocols for URL input include ftp, https and s3. S3 bucket access
can be configured using the URL syntax and environment variables supported by
samtools/htslib.

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
