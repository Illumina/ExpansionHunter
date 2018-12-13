# VCF files output by Expansion Hunter

Expansion Hunter generates a separate VCF record for each repeat with
information about repeat's location and genotype. The records for non hom-ref
repeats are demarcated by `<STRn>` symbolic alleles where `n` is the number of
repeat units that the corresponding allele spans.

The header of the VCF file contains a detailed description of each record.

## Example

The following VCF entry describes the state of *C9orf72* repeat in a sample with
name/barcode LP6005616-DNA_A03.

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  LP6005616-DNA_A03
chr9    27573526        .       C       <STR2>,<STR349> .       PASS    SVTYPE=STR;END=27573544;REF=3;RL=18;RU=GGCCCC;REPID=ALS GT:SO:CN:CI:AD_SP:AD_FL:AD_IR   1/2:SPANNING/INREPEAT:2/349:2-2/323-376:19/0:3/6:0/459
```

This line tells us that first allele spans 2 repeat units while the second
allele spans 349 repeat units. The repeat unit is GGCCCC (`RU` INFO field), so
the sequence of the first allele is GGCCCCGGCCCC and the sequence of the second
allele is GGCCCC x 349. The repeat spans three repeat units in the reference
(`REF` INFO field). The length of the short allele was estimated from spanning
reads (`SPANNING`) while the length of the expanded allele was estimated from
in-repeat reads (`INREPEAT`). The confidence interval for the size of the
expanded allele is (323,376). There are 19 spanning and 3 flanking reads
consistent with the repeat allele of size 2 (that is 19 reads fully contain the
repeat of size 2 and 2 flanking reads overlap at most 2 repeat units). Also,
there are 6 flanking and 459 in-repeat reads consistent with the repeat allele
of size 349.
