# Repeat catalogs

## Overview

A repeat catalog is an essential component of Expansion Hunter's input. It
specifies reference coordinates and structure of each repeat region that the
program will analyze. This document describes the format of files that store
repeat catalogs. Users that are considering creating new or modifying existing
catalog files should read this document carefully. Even a minor mistake in
specifying a repeat region could lead to a significant drop in genotyping
accuracy.

## Repeat catalog files

A repeat catalog file is a JSON array whose entries specify individual repeat
regions. Here is an example of a catalog consisting of three repeat regions:

```json
[
{
    "RepeatId": "DMPK",
    "RepeatUnit": "CAG",
    "ReferenceLocus": "19:46273463-46273522",
    "RepeatStatus": "common"
},
{
    "RepeatId": "FMR1",
    "RepeatUnit": "CGG",
    "ReferenceLocus": "X:146993569-146993628",
    "RepeatStatus": "rare",
    "OfftargetLoci": [
        "12:7781291-7781350",
        "12:125052155-125052156",
        "16:25703614-25703635",
        "16:28074517-28074518",
        "17:30814025-30814026",
        "17:64298468-64298469",
        "19:2015525-2015526",
        "2:87141541-87141618",
        "2:92230910-92230911",
        "2:211036021-211036032",
        "2:225449879-225449880",
        "20:30865501-30865516",
        "5:443335-443364",
        "7:20824940-20824941",
        "7:100271438-100271439",
        "7:104654598-104654599",
        "7:143059854-143059855",
        "9:100616696-100616697",
        "X:20009037-20009046"
    ]
},
{
    "RegionId": "HTT",
    "RegionStructure": "(CAG)CAACAG(CCG)",
    "RepeatIds": ["HTT_CAG", "HTT_CCG"],
    "ReferenceLoci": ["4:3076604-3076660", "4:3076667-3076693"],
    "RepeatStatuses": ["common", "common"]
}
]
```

The first entry in the above array specifies a repeat region consisting of a
single short tandem repeat. This repeat has identifier DMPK (field `RepeatId`)
and is comprised of repetitions of CAG repeat unit (field `RepeatUnit`). The
reference locus field specifies the exact coordinates of the repeat in the
reference (field `ReferenceLocus`). This repeat is "common" (field 
`RepeatStatus`) meaning that we expect the genome to contain multiple long 
repeats (whose size is close to fragment length and longer) with this repeat
unit. For "common" repeats Expansion Hunter limits the types of reads that are 
used used to infer the size of the repeat. As a result, common repeats are
genotyped up to a fragment length and if a repeat is reported to have a size
estimate close to the fragment length than this number should be treated as a
lower bound for its true size.

The second entry is similar to the first. The only conceptual difference is that
the status of the second repeat is set to "rare". This means that the user is 
*confident* that there are no other long (fragment length or longer) repeats 
with the same repeat unit elsewhere in the genome. For rare repeats, off-target
loci (field `OffTargetLoci`) specify regions of the genome that may contain
misaligned reads.

The final entry specifies a repeat region containing multiple short tandem 
repeats in close proximity to one another. As evidenced by the field names,
in such cases we make a distinction between regions containing the repeats and
the constituent repeats. The identifier of the entire region is HTT (field 
`RegionId`) while the identifiers of the two constituent repeats are HTT_CAG and
HTT_CCG (field `RepeatIds`). The structure of this region is encoded by a string
(CAG)CAACAG(CCG) where the short tandem repeats are specified by their repeat
unit enclosed in parenthesis with sequence between parentheses corresponding to
a non-repetitive interruption. Fields `ReferenceLoci` and `RepeatStatuses`
specify the exact reference coordinates and status of each repeat respectively 
(the entries have the same meaning as `ReferenceLocus` and `RepeatStatus` 
fields for regions containing a single repeat).

The following two sections give a detailed description of records specifying
repeats containing single and multiple short tandem repeats.

## Defining regions containing a single short tandem repeat

When a region contains a single short tandem repeat there is no difference
between the repeat and the region containing it. So field names refer only to
the repeat as shown in the table below.

| Field          | Description                                                                                |
|----------------|--------------------------------------------------------------------------------------------|
| RepeatId       | A unique identifier of the short tandem repeat                                             |
| RepeatUnit     | The repeat unit in the reference orientation                                               |
| ReferenceLocus | 1-based reference coordinates of the repeat (`chrom:start-end`)                            |
| RepeatStatus   | "common" if genome contains multiple long repeats with same repeat unit; "rare" otherwise  |
| OfftargetLoci  | Array of regions where informative reads may misalign (only for "rare" repeats)            |

## Defining regions containing multiple tandem repeats 

When a region contains multiple short tandem repeats, the corresponding record
must contain information about the region itself as well as about each
constituent repeat.

| Field           | Description                                                                    |
|-----------------|--------------------------------------------------------------------------------|
| RegionId        | A unique identifier of the repeat region                                       |
| RegionStructure | Region's structure encoded by a string*                                        |
| RepeatIds       | Array of unique identifiers for each constituent repeat                        |
| ReferenceLoci   | Array of reference coordinates for each repeat (see ReferenceLoci field above) |
| RepeatStatuses  | Array of repeat statuses for each repeat (see RepeatStatus field above)        |
| OfftargetLoci   | Array of additional regions with informative reads (only for "rare" repeats)   |

(*) The structure of this region is encoded by a string in the format
"(`repeat_unit`)`interruption`(`repeat_unit`) ... (`repeat_unit`)" that defines
multiple short tandem repeats and interrupting sequence between them (if any).
The short tandem repeats are defined by their repeat unit enclosed in
parenthesis with sequences between parentheses corresponding to interruptions.

## A note on creating custom repeat catalogs

Creating custom repeat catalogs is relatively straightforward for "common"
repeats. Defining "rare" repeats is much harder because some data analysis is
required to prove the the repeat is indeed rare. Users who are looking to define
custom repeat catalogs are encouraged to contact the developers for assistance.

It is important to specify `ReferenceLocus` in such a way that it encompasses
the full reference repeat sequence. That is, the sequence corresponding to
`ReferenceLocus` should (a) start and end with a perfect match to the repeat
unit and (b) the sequences adjacent to the repeat on both sides should be
distinct from the repeat unit.

One can use samtools to confirm that `ReferenceLocus` is specified correctly.
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
