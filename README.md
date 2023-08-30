# stantlib :toolbox: :microbe: :dna:
_Tom's toolbox for microbial sequence analysis_

 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8300404.svg)](https://doi.org/10.5281/zenodo.8300404)

## Installation :wrench:
Currently, there is no installation script. To use the toolbox, you need to clone the repository and add the path to 
your `PYTHONPATH` environment variable. For example, if you cloned the repository to `~/stantlib`, you would add the 
following line to your `~/.bashrc` file:

```bash
export PYTHONPATH=$PYTHONPATH:~/stantlib
```

## Overview :mag:
This repository contains a collection of scripts and tools I use for microbial sequence analysis. Where possible
I've written them in pure-Python which makes them fast and easy to use, but I've had to use dependencies for some
such as BioPython for parsing GenBank files :sweat:.

These scripts share common functions and classes which will eventually become a package, but for now they are just
standalone scripts. They are designed to be easy to use, and most accept standard input and pipe output which allows
them to be chained together.

Here is a brief overview of the main scripts:
- `fsk.py` - **F**asta **S**wiss-army **K**nife, a pure-Python tool to flatten, filter, and manipulate FASTA sequences 
from stdin. Regex can be used to filter sequences by accession/header. There is also simple utility to translate DNA 
sequences.
- `genbank2seq.py` - Extracts sequences from a GenBank file and writes them to stdout in FASTA format. Can be used to 
extract whole loci or features. If features are extracted, there is an option to extract the translation, either using
the `translation` qualifier or by translating the nucleotide sequence on-the-fly.
- `genbank_slicer.py` - Selects parts of a GenBank file and writes it to stdout in GenBank or FASTA format. Selection
methods include `--slice` with multiple `start:end` coordinates; `--feature` to select features with a regex or
`--record` to select records with a regex. The regex selection is performed on the string representation of the
feature/record, so it is very flexible!
- `gene_hotspots.py` - Looks for hotspots of metric windows across GenBank records and reports overlapping features.
This is useful for finding genes which are in areas of abnormal GC content, or which are in regions of high entropy.
- `gene_breaks.py` - Aligns an assembly to the reference and reports the positions of genes which are split across
contigs. This is useful for finding genes which are split across contigs, which can be a sign of misassembly.
- `gc_filter.py` - Calculate the GC content distribution across FASTQ reads and filter them.
- `clean_pairs.py` - Takes two paired FASTQ files and removes reads which are unpaired.
- `biopigz.py` - G(un)zip multiple files in parallel using `multiprocessing`. This is useful for unzipping many
small/medium genome assembly fasta files, in which case it is seems to be faster than `pigz`.
- `xccessions.py` - Create a regex from fixed-length accession numbers. This is useful for filtering sequences by
a header regex.

## Dependencies :package:
- Most of these tools require Python >= 3.8, and generally any Python version below this is not supported.
- Some scripts require `BioPython` which can be installed with `pip install biopython` or `conda install -c conda-forge
biopython`. 
- Some scripts also require external programs such as `minimap2`, but they will check for this and warn you if anything
is not installed.
