# Overview
Here, we describe the RNA-RNA Interaction Scan (rriScan), an universal software for analyzing RNA-RNA interactome sequencing data, such as `MARIO`, `PARIS`, `PARIS2` and `LIGR-seq` data.

- [System requirements](#System-Requirements)
- [Run time](#Run-time)
- [Installation](#Installation)
- [Input](#Input)
- [Output](#Output)
- [Basic usage](#Basic-usage)
- [Run rriScan on testing dataset](#Run-rriScan-on-testing-dataset)
- [Acknowledgements](#Acknowledgements)
- [Contact](#Contact)

# System requirements
The software package was tested on Linux system with RAM 15GB and CPU:20+ cores.

# Run time
Our software were written in C and C++ codes, so it has a supuer efficency when analyzing the RNA-RNA interactome sequencing data. After testing, most of the tasks can be finished within 10 minutes.

# Installation
* It's very easy to install rriScan on a linux server with following commands:
    ```bash
    # suppose your path to install software is /username/software
    
    cd /username/software
    
    git clone https://github.com/kerenzhou062/rriScan.git
    
    cd ./rriScan
    
    sh install.sh
    
    export PATH=$PATH:/username/software/rriScan/bin
    ```

# Input
* A genome sequences in fasta format is required by `--fa` argument, which chromsome's name must be the same with input `bam` file.

* An fai index file is required by `--fai` argument, which can be generated by using `samtools`:

    ```bash
    # suppose your input genome fasta file is hg38.fa
    
    samtools faidx hg38.fa
    ```

* A `bam` file contains sequence alignment data is required by `--bam` argument, which can be generated by `STAR` software and stored as `Aligned.sortedByCoord.out.bam` file.

* A `juction` file contains junctions information, which can be generated by `STAR` software with `--chimOutType Junctions` parameter and stored as `Chimeric.out.junction` file.

* A `fastq` or `fasta` file contains sequencing reads is optional.

# Output
Here's the description of columns in the outputs:

| Column name          | Description
| -----------          |----------
| `lChrom`             | chromosome name of left pair
| `lChromStart`        | start coordinate of left pair (0-base)
| `lChromEnd`          | end coordinate of left pair
| `lName`              | name of left pair
| `lScore`             | score of left pair
| `lStrand`            | strand of left pair
| `rChrom`             | chromosome name of right pair
| `rChromStart`        | start coordinate of right pair (0-base)
| `rChromEnd`          | end coordinate of right pair
| `rName`              | name of right pair
| `rScore`             | score of right pair
| `rStrand`            | strand of right pair
| `lociNum`            | chromosome name
| `gapDist`            | gap distance between pairs if they are in the same chromosome
| `readSeq`            | sequencing reads
| `chimericSeq`        | full sequence of the chimeric
| `chimericStruct`     | predict structure of chimeric
| `MFE`                | minimum free energy
| `rriType`            | type of RNA-RNA interaction
| `lAlignSeq`          | aligned sequence of left pair
| `pairs`              | base pairings
| `rAlignSeq`          | aligned sequence of right pair
| `pairNum`            | the maximum continuous perfect pairings
| `alignScore`         | Smith-Waterman score
| `loReadNum`          | read number of left pair
| `roReadNum`          | read number of right pair

# Basic Usage
The available options of rriScan are as follow:

```shell
Usage:  rriScan [options] --fa <fasta file> --fai <fai file> --bam <mapped alignments> --jun <junctions>
File format for mapped alignments is BAM
[options]
-v/--verbose                   : verbose information
-V/--version                   : rriScan version
-h/--help                      : help informations
-S/--small                     : small genome
--fa                           : genome FASTA file. [required]
--fai                          : genome fai file, an index file for fasta file. [required]
--bam                          : alignment file, BAM format. [required]
--jun                          : junction file from STAR software, junction format. [required]
--read                         : read file[fastq or fasta]. [optional]
-o/--output <string>           : output file
-l/--min-seg-len <int>         : minimum length of segments in a chimera [default>=15]
-m/--min-read-number <int>     : minimum read number for chimera [default>=1]
-M/--max-mfe <double>          : maximum MFE in duplex[default<=-5.0]
-p/--min-pair <int>            : minimum pair number in duplex [default>=0]
-s/--min-score <int>           : minimum alignment score in duplex [default>=5]
-g/--min-gap <int>             : minimum gaps between two segments [default>=1]
```

# Run rriScan on testing dataset
Please [check here](test_data/README.md) to learn how to run rriScan on a testing dataset.

# Acknowledgements
Thanks a lot to everyone who contributed to the public codes and libraries (e.g. BamTools) used by rriScan.

# Contact
* Jian-Hua Yang <yangjh7@mail.sysu.edu.cn>, RNA Information Center, School of Life Sciences, Sun Yat-Sen University<BR>
* Keren Zhou <kzhou@coh.org>, Department of Systems Biology, Beckman Research Institute of City of Hope, Monrovia, CA, USA<BR>
