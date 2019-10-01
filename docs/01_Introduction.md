# A guide to Expansion Hunter software

Expansion Hunter is a tool for targeted genotyping of short tandem repeats and
flanking variants. It operates by performing a targeted search through a
BAM/CRAM file for reads that span, flank, and are fully contained in each
repeat. Newer versions of the program provide limited support for insertions,
deletions, and substitutions.

In order to use Expansion Hunter you need to (1) download the latest release or
build the program from source, (2) obtain a set of files specifying repeat
regions of interest, and (3) have a BAM or a CRAM file with alignments of reads
from a PCR-free WGS sample.

The following sections will help you to get started.

* [Installation](02_Installation.md)
* [Usage](03_Usage.md)
* [Input variant catalogs](04_VariantCatalogFiles.md)
* [Output JSON files](05_OutputJsonFiles.md)
* [Output VCF files](06_OutputVcfFiles.md)

Important limitations:

- Regions whose coverage depth is below 10x are skipped during the analysis.
  This is done because the STR genotyping accuracy drops significantly when
  the coverage dips below 10-15x.
- ExpansionHunter does not call novel variants; it only genotypes variants
  defined in the variant catalog file. In particular, the program does not
  report interruptions in the repeat sequence.
- The current release (v3.x.x) forgoes analysis of unaligned reads. Although
  some informative read pairs are occasionally unaligned, according to our
  benchmarks these reads are very unlikely to change pathogenic/normal
  classification of the currently-covered repeats. We encourage you to contact
  the authors if you are working with the regions for which this is a concern.
