# Expansion Hunter: a tool for estimating repeat sizes

There are a number of regions in the human genome consisting of repetitions of
short unit sequence (commonly a trimer). Such repeat regions can expand to a
size much larger than the read length and thereby cause a disease.
[Fragile X Syndrome](https://en.wikipedia.org/wiki/Fragile_X_syndrome),
[ALS](https://en.wikipedia.org/wiki/Amyotrophic_lateral_sclerosis), and
[Huntington's Disease](https://en.wikipedia.org/wiki/Huntington%27s_disease)
are well known examples.

Expansion Hunter aims to estimate sizes of such repeats by performing a targeted
search through a BAM/CRAM file for reads that span, flank, and are fully
contained in each repeat.

Linux and macOS operating systems are currently supported.


## License

Expansion Hunter is provided under the terms and conditions of the
[GPLv3 license](LICENSE.txt). It relies on several third party packages provided
under other open source licenses, please see [COPYRIGHT.txt](COPYRIGHT.txt) for
additional details.


## Documentation

Installation instructions, usage guide, and description of file formats are
contained in the [docs folder](docs/01_Introduction.md).


## Method

The detailed description of the method can be found here:

Dolzhenko and others, [Detection of long repeat expansions from PCR-free
whole-genome sequence data](http://genome.cshlp.org/content/27/11/1895), Genome
Research 2017
