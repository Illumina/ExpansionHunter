# Outputs

VCF output
----------

The VCF output file contains the following fields.

 Field  | Description
--------|----------------------------------------------------------------------------------
 CHROM  | Chromosome identifier
 POS    | Position of the first base before the repeat region in the reference
 ID     | Always `.`
 REF    | The reference base at position POS
 ALT    | List of repeat alleles in format `<STRn>` where n is the number of repeat units
 QUAL   | Always `.`
 FILTER | Always PASS

These are followed by the INFO fields.

Field    | Description
---------|---------------------------------------------------------
 SVTYPE  | Always STR
 END     | Position of the last base of the repeat region in the reference
 REF     | Number of repeat units spanned by the repeat in the reference
 RL      | Reference length in bp
 RU      | Repeat unit in the reference orientation
 REPID   | Repeat id from the repeat-specification file

And finally the sample fields.

 Field | Description
:-----:|---------------------------------------------------------
 GT    | Genotype
 SO    | Type of reads that support the allele; can be SPANNING, FLANKING, or INREPEAT meaning that the reads span, flank, or are fully contained in the repeat
 CN    | Allele copy number
 CI    | Confidence interval for CN
 AD_SP | Number of spanning reads consistent with the allele
 AD_FL | Number of flanking reads consistent with the allele
 AD_IR | Number of in-repeat reads consistent with the allele

For example, the following VCF entry describes the state of 
*C9orf72* repeat in a sample with ID LP6005616-DNA_A03.

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  LP6005616-DNA_A03
chr9    27573526        .       C       <STR2>,<STR349> .       PASS    SVTYPE=STR;END=27573544;REF=3;RL=18;RU=GGCCCC;REPID=ALS GT:SO:CN:CI:AD_SP:AD_FL:AD_IR   1/2:SPANNING/INREPEAT:2/349:2-2/323-376:19/0:3/6:0/459
```

This line tells us that first allele spans 2 repeat units while 
the second allele spans 349 repeat units. The repeat unit is 
GGCCCC (`RU` INFO field), so the sequence of the first allele 
is GGCCCCGGCCCC and the sequence of the second allele is GGCCCC x 349. 
The repeat spans three repeat units in the reference (`REF` INFO 
field). The length of the short allele was estimated from spanning 
reads (`SPANNING`) while the length of the expanded allele was 
estimated from in-repeat reads (`INREPEAT`). The confidence interval 
for the size of the expanded allele is (323,376). There are 19 
spanning and 3 flanking reads consistent with the repeat allele of 
size 2 (that is 19 reads fully contain the repeat of size 2 and 2 
flanking reads overlap at most 2 repeat units). Also, there are 6 
flanking and 459 in-repeat reads consistent with the repeat allele 
of size 349.

JSON output
-----------

Each repeat region is represented by a JSON object. Here is the list 
of fields contained in each object and their description.

 Field             | Description
-------------------|---------------------------------------------------------
 RepeatId          | A unique string identifier of the repeat region
 TargetRegion      | 1-based coordinates of the repeat region in the reference, specified by `chrom:start-end`
 RepeatUnit        | The unit of the repeat
 Genotype          | Repeat genotype; a pair of repeat sizes separated by `/` for diploid chromosomes and a single repeat size for haploid chromosomes
 GenotypeCi        | Confidence interval for the size of each repeat allele
 GenotypeSupport   | The number of spanning, flanking, and in-repeat reads (in this order) consistent with each repeat allele
 IrrCount          | The total number of identified in-repeat reads
 AnchoredIrrCount  | The number of in-repeat reads anchored by their mates to the repeat region
 OffTargetRegionIrrCounts | Is an object consisting of region/count pairs with each count giving the number of in-repeat reads found in the corresponding region
 UnalignedIrrCount | The number in-repeat reads found in the unaligned section of the BAM/CRAM file
 RepeatSizes       | Describes the type of the reads (`Source`), their number (`NumSupportingReads`), and how many repeat units they contain (`Size`); for flanking reads, `NumSupportingReads` gives the size of the longest flanking read identified

The log file
------------

The log file contains alignments of spanning and flanking reads and sequences of in-repeat reads.

The top-level keys are repeat region identifiers and the values correspond to repeat alleles. Each allele is described by `<READ_TYPE>_<SIZE>` where `READ_TYPE` can be `SPANNING`, `FLANKING`, or `INREPEAT` and `SIZE` is the size of the repeat in repeat units. For each allele estimated from spanning/flanking reads, the log file lists the name and the alignment of each spanning/flanking read respectively. The low quality bases are printed in lowercase. For `INREPEAT` alleles, the log file lists the sequences of in-repeat reads and, if appropriate, anchors. One of the records is labeled as `FLANKING`; it contains alignments of flanking reads that likely came from one of the `SPANNING` alleles.