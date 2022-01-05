# rnaseq

Proof of concept of a RNA-Seq pipeline from .fastq reads to count matrix (including quality control) with Nextflow.

## Prerequisites

* Unix-like OS (Linux, macOS, etc.)
* [Java](https://openjdk.java.net) version 8
* [Docker](https://docs.docker.com/engine/install/) engine 1.10.x (or later)

### Additional necessary files

* Prebuild Hisat2 index for H. sapiens, release GRCh38

```
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_snptran.tar.gz
```

* Gencode GTF file, release 38 (GRCh38.p13)

```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz
```

* USCS BED file, assembly GRCh38/hg38 track GENCODE V38

```
http://genome.ucsc.edu/cgi-bin/hgTables
```

## Table of Contents

* [Quick start](#Quick-start)
* [Installation](#Installation)
* [Arguments](#Arguments)
* [Documentation](#Documentation)

### Quick start

### Installation

### Arguments

### Documentation

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:
