# Variant catalogs

## Overview

A variant catalog is an essential component of Expansion Hunter's input. It
specifies reference coordinates and structure of each locus that the program
will analyze. Although the loci usually correspond to repeats, they could also
contain other classes of variants such as insertions, deletions, and sequence
swaps.

This document describes the format of files that store repeat catalogs. Users
that are considering creating new or modifying existing catalog files should
read this document carefully. Even a minor mistake in locus definition could
lead to a significant drop in genotyping accuracy.


## Variant catalog files

A variant catalog file is a JSON array whose entries specify individual loci
that the program will analyze. Here is an example of a catalog consisting of
three loci containing repeats.

```json
[
{
    "LocusId": "DMPK",
    "LocusStructure": "(CAG)*",
    "ReferenceRegion": "19:46273462-46273522",
    "VariantType": "Repeat"
},
{
    "LocusId": "FMR1",
    "LocusStructure": "(CGG)*",
    "ReferenceRegion": "X:146993568-146993628",
    "VariantType": "RareRepeat",
    "OfftargetRegions": [
        "12:7781290-7781350",
        "12:125052154-125052156",
        "16:25703613-25703635",
        "16:28074516-28074518",
        "17:30814024-30814026",
        "17:64298467-64298469",
        "19:2015524-2015526",
        "2:87141540-87141618",
        "2:92230909-92230911",
        "2:211036020-211036032",
        "2:225449878-225449880",
        "20:30865500-30865516",
        "5:443334-443364",
        "7:20824939-20824941",
        "7:100271437-100271439",
        "7:104654597-104654599",
        "7:143059853-143059855",
        "9:100616695-100616697",
        "X:20009036-20009046"
    ]
},
{
    "LocusId": "HTT",
    "LocusStructure": "(CAG)*CAACAG(CCG)*",
    "ReferenceRegion": ["4:3076604-3076660", "4:3076666-3076693"],
    "VariantType": ["Repeat", "Repeat"]
}
]
```

The first entry in this catalog specifies a locus containing a single short
tandem repeat. The identifier of this locus is DMPK (field `LocusId`). The
regular expression `(CAG)*` means that it is comprised of zero or more
repetitions of the CAG repeat unit (field `LocusStructure`). The reference
coordinates of this repeat are 19:46273462-46273522 (field  `ReferenceRegion`).
The `VariantType` field specifies that it is an ordinary STR meaning that we
expect the genome to contain multiple long repeats (whose size is close to
fragment length and longer) with this repeat unit. For ordinary `Repeat`s
Expansion Hunter limits the types of reads that are used used to infer the size
of the repeat. As a result, regular repeats are genotyped up to the fragment
length and if a repeat is reported to have a size estimate close to the fragment
length then this number should be treated as a lower bound for its true size.

The second entry is similar to the first. The only difference is that the
variant type of the second repeat is set to `RareRepeat`. This means that the
user is *confident* that there are no other long (fragment length or longer)
repeats with the same repeat unit elsewhere in the genome. This information
permits the program to use additional read-level evidence to potentially
estimate the length of the repeat past the fragment length. For "rare" repeats,
off-target regions (field `OfftargetRegions`) specify regions of the genome that
may contain misaligned reads.

The final entry describes a repeat region containing multiple short tandem
repeats in close proximity to each other. The regular expression
`(CAG)*CAACAG(CCG)*` specifies that this region consists of two short tandem
repeats with repeat units CAG and CCG separated by the sequence CAACAG.
Fields `ReferenceRegion` and `VariantStatus` contain reference region and status
of each constituent repeat. By default, the program assigns an identifier to
each variant consisting of the locus id and reference region. So the two repeats
receive ids HTT_4:3076603-3076660 and HTT_4:3076666-3076693 respectively. An
optional field `VariantId` allows to assign custom variant ids to each
variant.

The following section describes loci-specification records that the catalogs are
comprised of.


## Structure of a locus-specification record

When locus contains a single variant, there is no difference between the
variant and the locus containing it. So field names refer to the variant
itself.

* `LocusId` Unique identifier of the entire locus

* `LocusStructure` Regular expression defining the structure of the locus. When
  the locus contains multiple variants, ReferenceRegion and VariantStatus are
  arrays with associated information for each variant in the same order.

* `ReferenceRegion` 0-based half open reference coordinates of the variant
  formatted as `chrom:start-end`.

* `VariantType` Can be either `Repeat`, `RareRepeat`, or `SmallVariant`
  with the latter corresponding to insertions deletions or sequence swaps.

* `VariantId` Optional array of unique variant ids. If missing, variant ids
  are synthesized according to this rule: If there is only one variant in
  a locus then it gets the same id as the locus itself. If locus contains
  multiple variants, each one of them gets id of the form `<LocusId>_<ReferenceRegionOfTheVariant>`.

* `OfftargetRegions` Array of regions where informative reads may misalign;
   only used for variants of type `RareRepeat`.


## Using regular expressions to define locus structure

ExpansionHunter supports a very limited subset of regular expressions to
define the structure of each locus. These expressions can consist of
sub-expressions listed in the table bellow, possibly separated by
interrupting DNA sequences.


| Variant                                            | Regular expression |
|----------------------------------------------------|--------------------|
| Short tandem repeat that can occur 0 or more times | (CCG)*             |
| Short tandem repeat that can occur 1 or more times | (CCG)+             |
| Single nucleotide variant                          | (C\|T)             |
| Sequence swap                                      | (CAGT\|CGTTG)      |
| Deletion or insertion	                             | (CTGGC)\?          |


For example, a CAG repeat flanked by a CAG/CAT swap	is defined by expression
(CAG)+CTGT(CAG|CAT).


## A note on creating custom variant catalogs

Creating custom variant catalogs is relatively straightforward for "common"
variants. Defining "rare" variants is much harder because some data analysis is
required to prove that the variant is indeed rare. Users who are looking to
define custom catalogs are encouraged to contact the developers for assistance.
