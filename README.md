# rnaseq

Proof of concept of a RNA-Seq pipeline from reads to count matrix (including quality control) with Nextflow and subsequent RNA-Seq analysis in R.

## Prerequisites

* Unix-like OS (Linux, macOS, etc.)
* [Java](https://openjdk.java.net) version 8
* [Docker](https://docs.docker.com/engine/install/) engine 1.10.x (or later)


### Necessary files

* Reads to be mapped must be stored in compressed `.fastq.gz` file format in folder `data`

### Additional necessary files

These additional 3 files must be stored in folder `data`:

* Prebuild Hisat2 index for H. sapiens, release GRCh38

```
https://genome-idx.s3.amazonaws.com/hisat/grch38_snptran.tar.gz
```

* Gencode GTF file, release 38 (GRCh38.p13)

```
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz
```

* USCS BED file, assembly GRCh38/hg38, track GENCODE V38

```
http://genome.ucsc.edu/cgi-bin/hgTables
```
The BED file must be stored in `*.annotation.bed.gz` file format.

## Table of Contents

* [Quick start](#Quick-start)
* [Installation](#Installation)
* [Arguments](#Arguments)
* [Documentation](#Documentation)

## Quick start

Example run:
```
nextflow run main.nf
```

The above example uses default parameter `params.reads` for single-end reads:
```
nextflow run main.nf --reads "data/*.fastq.gz"
```

For paired-end reads, additionally parameter `params.singleEnd` in `nextflow.config` must be changed to `false`. Then the input command must be:

```
nextflow run main.nf --reads "data/*_{1,2}*.fastq.gz"
```

Optionally, you can specify the Nextflow output directory with flag `--outdir <folder>`. By default, all resulting files will be saved in folder `output` and folder `info` will contain all information about the last run nextflow session.

## Installation

Clone this repository with the following command:

```
git clone https://github.com/maxgreil/rnaseq && cd rnaseq
```

Then, install Nextflow by using the following command:

```
curl https://get.nextflow.io | bash
```

The above snippet creates the `nextflow` launcher in the current directory.

Finally pull the following Docker container:

```
docker pull maxgreil/rnaseq
```

Alternatively, you can build the Docker Image yourself using the following command:

```
cd docker && docker image build . -t maxgreil/rnaseq
```

## Arguments

### Optional Arguments

| Argument  | Usage                            | Description                                                          |
|-----------|----------------------------------|----------------------------------------------------------------------|
| --reads| \<files\>                           | Directory and glob pattern of input files|
| --outdir  | \<folder\>                       | Directory to save output files                                    |

## Documentation

This pipeline is designed to:

* map given reads to a genome
* create a count matrix of mapped reads for subsequent RNA-Seq analysis
* do a quality control of the created files

### Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

1. [hisat2](http://daehwankimlab.github.io/hisat2/) - map given reads to genome
2. [samtools](http://www.htslib.org/) - create sorted BAM files from HISAT2 SAM files
3. [featureCounts](http://subread.sourceforge.net/) - count mapped reads to genomic features (exons)
4. [preseq](http://smithlabresearch.org/software/preseq/) -  predict and estimate the complexity of genomic sequencing library
5. [reseqc](http://rseqc.sourceforge.net/) - comprehensive evaluation of used RNA-Seq data
6. [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - BAM quality control
7. [MultiQC](https://multiqc.info) - aggregate report, describing results of the whole pipeline
