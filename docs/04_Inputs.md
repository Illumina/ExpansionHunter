# Inputs

Repeat-specification files
--------------------------

Repeat-specification files define repeat regions for Expansion Hunter to
analyze. Each repeat-specification file is a JSON file containing a single 
object that consists of name/value pairs. Sample repeat-specification for 
some pathogenic repeats are contained in `data/repeat-specs/` directory.

| Field            | Description                                                                               |
|------------------|-------------------------------------------------------------------------------------------|
| RepeatId         | A unique string identifier of the repeat region.                                          |
| RepeatUnit       | The repeat unit (in reference orientation) that the repeat is comprised of.               |
| CommonUnit       | If true, only anchored IRRs are used to estimate the long repeat sizes* (Default: false). |
| TargetRegion     | 1-based coordinates of the repeat region in the reference, specified as `chrom:start-end`.|
| OffTargetRegions | An array of off-target regions. This field is optional and used only if `CommonUnit` is set to `false`.                                                         |

\* Should be set to true if repeats longer than the read length and having 
the same repeat unit are expected to occur elsewhere in the genome.

A note on creating custom repeat-specification files
----------------------------------------------------

Creating specification files for new repeat regions is easy: Just use one of 
the provided specification files as a template.

It is important to specify `TargetRegion` in such a way that it encompasses 
the full reference repeat sequence. That is, the reference sequence of the 
`TargetRegion` should (a) start and end with a perfect match to the repeat 
unit and (b) the sequences adjacent to the repeat on both sides should be 
distinct from the repeat unit.    

One can use `samtools` to confirm that `TargetRegion` is specified correctly. 
For example, the target region of *C9orf72* repeat (chr9:27573527-27573544 in 
hg19) is a perfect repetition of GGCCCC hexamer:
```bash
  $ samtools faidx hg19_ref.fa chr9:27573527-27573544
    >chr9:27573527-27573544
    ggccccggccccggcccc
```
while the sequences adjacent to the left and right sides of the repeat do not 
end and start (respectively) with the repeat unit hexamer:
```bash
  $ samtools faidx hg19_ref.fa chr9:27573507-27573526
    >chr9:27573507-27573526
    gcccgccccgaccacgcccc

  $ samtools faidx hg19_ref.fa chr9:27573545-27573564
    >chr9:27573545-27573564
    TAGCGCGCGACTCCTGAGTT
```
