# Expansion Hunter: a tool for estimating repeat sizes


## Overview of long-term OSS maintenance fork

ExpansionHunter was converted from open (Apache 2 or GPL 3) to Polyform Strict
licensing after v5.0.0. This fork is designed to preserve the open source
codebase used in the releases up to and including v5.0.0, and to provide
a repo for maintenance patches moving forward.


### Licensing

As of v5.0.0, licensing of the codebase was ambiguous. The LICENSE.txt
file included in the repository was Apache 2, and all source contributions
outside of third party content are Apache 2. However, the COPYRIGHT.txt file
in the repo lists the package as being GPL 3. The LICENSE and COPYRIGHT
files are thus in conflict in v5.0.0. Based on file editing times and the
preponderance of code contributions, it appears in this case that the GPL
addition was made in error, and in fact the code base was predominantly
made with the intention of creating an Apache 2 project. As such, to try
to create an internally consistent license state that respects the history
of the repository as best as possible, the COPYRIGHT text is being adjusted
to be consistent with the existing LICENSE file for all releases after v5.0.0.


## Original overview from ExpansionHunter repo



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
[Apache License Version 2.0](LICENSE.txt). It relies on several third party
packages provided under other open source licenses, please see
[COPYRIGHT.txt](COPYRIGHT.txt) for additional details.


## Documentation

Installation instructions, usage guide, and description of file formats are
contained in the [docs folder](docs/01_Introduction.md).


## Method

The method is described in the following papers:

- Egor Dolzhenko, Joke van Vugt, Richard Shaw, Mitch Bekritsky, and others,
  [Detection of long repeat expansions from PCR-free whole-genome sequence data](http://genome.cshlp.org/content/27/11/1895),
  Genome Research 2017

- Egor Dolzhenko, Viraj Deshpande, Felix Schlesinger, Peter Krusche, Roman Petrovski, and others,
[ExpansionHunter: A sequence-graph based tool to analyze variation in short tandem repeat regions](https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btz431/5499079),
Bioinformatics 2019
