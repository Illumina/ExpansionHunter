
Expansion Hunter can be run on the example BAMlet using the following command
after replacing "hg19_ref.fa" with a path to hg19 reference FASTA file and 
adjusting paths in the remaining command line arguments as approprite.

Note that "read-depth" must be specified when running Expansion Hunter on an
incomplete BAM file.

ExpansionHunter \
  --bam bamlets/bamlet.bam \
  --ref-fasta hg19_ref.fa \
  --repeat-specs ../../repeat-specs/hg19/ \
  --vcf output/bamlet.vcf \
  --json output/bamlet.json \
  --log output/bamlet.log \
  --read-depth 30
