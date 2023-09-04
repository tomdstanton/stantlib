# stantlib :toolbox: :microbe: :dna:
_Tom's toolbox for microbial sequence analysis_

 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8300404.svg)](https://doi.org/10.5281/zenodo.8300404)

This repo contains a collection of scripts and tools I use for microbial sequence analysis.
They share common functions and classes which will eventually 
become a package, but for now they are just standalone scripts. **Please cite if you found anything useful** :pray:.

Please feel free to contribute, especially if you think you can make anything faster or more efficient :rocket:.
I may switch to using corresponding `gff + fasta` inputs instead of `GenBank` files if there is a demand to completely
remove `BioPython` dependencies.

## Installation :wrench:
Currently, there is no installation script. To use the toolbox, you need to clone the repository and add the path to 
your `PYTHONPATH` environment variable. For example, if you cloned the repository to `~/stantlib`, you would add the 
following line to your `~/.bashrc` file:

```bash
export PYTHONPATH=$PYTHONPATH:~/stantlib
```

## Dependencies :package:
- Most of these tools require Python >= 3.8 due to the use of the `walrus` operator `:=`, f-strings, and type hints.
- Some scripts require `BioPython` which can be installed with `pip install biopython` or `conda install -c conda-forge
biopython`. 
- Some scripts also require external programs such as `minimap2`, but they will check for this and warn you if anything
is not installed.


## Overview :mag:
Detailed usage of each script can be found by running it with the `-h` / `--help` flag; here is a brief description in 
no particular order:

### fsk :hammer_and_pick:
- **F**asta **S**wiss-army **K**nife to flatten, filter, and manipulate FASTA sequences from stdin.
- Regex can be used to filter sequences by accession/header. 
- There is also simple utility to translate DNA sequences.
- No dependencies!

```shell
# Extract ompC gene sequences from a collection of annotated genomes and translate them
$ cat bakta/*.ffn | fsk.py - --regex '(?i)ompc' -t -T 11 > OmpC.faa  # Note the python regex flag (?i) for case-insensitive
```

### gluggins :wine_glass:
- Parse Gubbins output as recombinant blocks of features using the reference sequence Genbank file.
- A block is defined as a region of the reference sequence with overlapping recombination events predicted by Gubbins.
- Blocks will only be constructed if they meet the selected critera and features can be reported by type (e.g. CDS, tRNA, rRNA, etc.).
- All features overlapping the block are reported, along with the recombination events which overlap them.
- Requires `BioPython`

```shell
# Explore which genes are in recombinant blocks
$ gluggins.py reference.gbff recombination_predictions.gff | head
reference	block_id	block_start	block_end	feature_type	feature_start	feature_end	feature_strand	feature_name	feature_locus_tag	feature_description	events
chr	Block_1	13311	20983	CDS	13311	14603	+	uhpC	CKEIMB_00065	MFS transporter	13765 17554 Node_46->Node_45 14 1228224,6293648
chr	Block_1	13311	20983	CDS	14655	15875	+		CKEIMB_00070	DUF3748 domain-containing protein	13765 17554 Node_46->Node_45 14 1228224,6293648
chr	Block_1	13311	20983	CDS	15877	16209	-	yceK	CKEIMB_00075	YceK/YidQ family lipoprotein	13765 17554 Node_46->Node_45 14 1228224,6293648
chr	Block_1	13311	20983	CDS	16524	16937	+	ibpA	CKEIMB_00080	heat shock chaperone IbpA	13765 17554 Node_46->Node_45 14 1228224,6293648
chr	Block_1	13311	20983	CDS	17054	17482	+	ibpB	CKEIMB_00085	heat shock chaperone IbpB	13765 17554 Node_46->Node_45 14 1228224,6293648
chr	Block_1	13311	20983	CDS	17619	18443	+	dinD	CKEIMB_00090	DNA damage-inducible protein D	18011 20813 Node_46->Node_45 34 1228224,6293648
```

### genbank_slicer :book::hocho:
- Selects parts of a GenBank file and writes it to stdout in GenBank or FASTA format. 
- Selection methods include `--slice` with multiple `start:end` coordinates; `--feature` to select features with a regex or
`--record` to select records with a regex. 
- The regex selection is performed on the string representation of the
feature/record (including the sequence/translation), so try and make your regex specific!
- Requires `BioPython`

```shell
# Extract the cps locus from a collection of annotated genomes using the feature regex
$ cat bakta/*.gbff | genbank_slicer.py - --feature 'galF|ugd' --merge 30000 > cps.gbk
```

### genbank2seq :book::dna:
- Extracts sequences from a GenBank file and writes them to stdout in FASTA format. 
- Can be used to extract whole loci or features. 
- If features are extracted, there is an option to extract the translation, either using
the `translation` qualifier or by translating the nucleotide sequence on-the-fly.
- Requires `BioPython`

```shell
# Extract the translation of all CDS features from a GenBank file and write them to a FASTA file
$ cat *.gbk | genbank2seq.py - --record feature --feature-type CDS --protein > hello.faa
```

### biopigz :pig:
- G(un)zip multiple files in parallel using `multiprocessing`. 
- This is useful for unzipping many small/medium genome assembly fasta files, in which case it is seems to be faster 
than `pigz` on my M1 chip.
- Processes files in place, adding/removing the '.gz' suffix.
- If stdin and not stdout, processed output will be written to [uuid.file(.gz)]
- Multiple files will be processed simultaneously using `--threads`
- No dependencies!

```shell
# Decompress all gzipped FASTA files in the current directory using 8 threads and keep the original files
$ biopigz.py *.fasta.gz --threads 8 --decompress --keep
```

### gc_filter :chart_with_upwards_trend::goal_net:
- Calculate the GC content distribution across FASTQ reads and filter them.
- No dependencies!

```shell
# Filter reads with a GC content that is 2 standard deviations from the mean
$ gc_filter.py reads.fastq --zscore 2 > filtered.fastq
# Use pigz to speed up (de)compression (remember biopigz is for lots of smallish files)
$ pigz -dc reads.fastq.gz | gc_filter.py - | pigz > filtered.fastq.gz
# Select your own threshold
$ gc_filter.py reads.fastq --upper 0.6 --lower 0.4 > filtered.fastq
```

### clean_pairs :sponge::couple:
- Takes two paired FASTQ files and removes reads which are unpaired.
- No dependencies!

```shell
$ clean_pairs.py reads_{1,2}.fastq.gz
```

### gene_gc :chart_with_upwards_trend::dna:
- Calculates the GC content of features in a GenBank file.
- Requires `BioPython`

```shell
# Calculate the GC content of all ncRNA features in a GenBank file
$ gene_gc.py genome.gbff -f ncRNA | head
reference	reference_gc	feature_type	feature_start	feature_end	feature_name	feature_locus_tag	feature_product	gc	gc_diff
chr	0.575	ncRNA	AsdA	CKEIMB_00010	antisense RNA of dnaA mRNA	510	591	0.617	0.042
chr	0.575	ncRNA	istR	CKEIMB_00165	istR Hfq binding RNA	32653	32784	0.489	-0.087
chr	0.575	ncRNA	STnc370	CKEIMB_00955	Enterobacterial sRNA STnc370	198260	198330	0.314	-0.261
chr	0.575	ncRNA	sok	CKEIMB_00960	sok antitoxin (CssrC)	198427	198580	0.451	-0.124
chr	0.575	ncRNA	RyhB	CKEIMB_01560	RyhB RNA	336435	336501	0.455	-0.121
```

### gene_hotspots :dna::fire:
- Looks for hotspots of metric windows across GenBank records and reports overlapping features.
- This is useful for finding genes which are in areas of abnormal GC content, or which are in regions of high entropy.
- Requires `BioPython`

```shell
 # Find genes which are in areas of entropy 3 or more standard deviations from the mean
 $ gene_hotspots.py genome.gbff --analysis entropy --zscore 3 | head
 reference	start	stop	analysis	value	feature_type	feature_start	feature_end	feature_name	feature_locus_tag	feature_description
chr	227250	227750	entropy_hotspot	1.89	CDS	227066	227969	dppC	CKEIMB_01100	dipeptide ABC transporter permease DppC
chr	237750	238250	entropy_hotspot	1.88	CDS	236370	240423	hemYx	CKEIMB_01135	cellulose synthase
chr	253750	254250	entropy_hotspot	1.89	CDS	252611	256091	bcsC	CKEIMB_01190	cellulose synthase complex outer membrane protein BcsC
chr	308000	308500	entropy_hotspot	1.88	CDS	306781	308992	zntA	CKEIMB_01415	Zn(II)/Cd(II)/Pb(II) translocating P-type ATPase ZntA
chr	389750	390250	entropy_hotspot	1.88	CDS	389565	390099	pilN	CKEIMB_01770	pilus assembly protein PilN
chr	389750	390250	entropy_hotspot	1.88	CDS	390088	390517		CKEIMB_01775	pilus assembly protein PilO
```
### gene_breaks :dna::broken_heart:
- Aligns an assembly to the reference and reports features which are split across contigs or missing.
- Requires `BioPython` and `minimap2`

```shell
# Compare a genome assembled from Illumina reads to a reference genome and report incomplete/missing CDS
$ gene_breaks.py completed_genome.gbff illumina_assembly.fasta --feature_type CDS | head
reference	reference_start	reference_end	contig	contig_start	contig_end	feature_type	feature_start	feature_end	feature_name	feature_locus_tag	feature_description
contig_1	1030282	1231824	5	8	201426	CDS	1029194	1030379	tuf	OAFDHN_05205	elongation factor Tu
contig_1	4958761	5153929	6	112	195179	CDS	4958460	4960620		OAFDHN_23805	Polysaccharide biosynthesis tyrosine autokinase
contig_1	4231604	4413295	7	3	181694	CDS	4231426	4232371	araC	OAFDHN_20265	Rhamnose utilization AraC family transcriptional regulator
contig_1	696620	863324	8	0	166704	CDS	696289	696883		OAFDHN_03405	Lipoprotein
contig_1	257038	415548	10	53	158563	CDS	415288	415987		OAFDHN_02010	ABC transporter permease
contig_1	4481764	4616495	12	346	135077	CDS	4481517	4482063		OAFDHN_21530	Thiolase-N domain-containing protein
contig_1	863468	998040	13	1	134457	CDS	863432	863897		OAFDHN_04270	OB-fold putative lipoprotein
```

### xccessions :name_badge:
- Create a regex from fixed-length accession numbers. 
- This is useful for filtering sequences by a header regex.
- No dependencies!

```shell
# Create a regex from a list of NCBI accessions
$ xccessions.py WP_000000000.1 WP_000000001.1 WP_000000002.1
WP_00000000[120].1
# Create a regex from a list of UniProt accessions
$ xccessions.py $(cat uniprot_accessions.txt)
```
### get_hmms :face_in_clouds:
- Download Hidden-Markov/Covariance models from public databases with accessions.
- Uses `multiprocessing` to download multiple models simultaneously.
- Requires `requests`

```shell
# Download 3 Rfam covariance models using 4 threads and use cm as the file extension
$ get_hmms.py RF00444 RF00445 RF00446 --extention cm --threads 4
```
