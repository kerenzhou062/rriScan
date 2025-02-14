# run rriScan on testing CLIP datasets
Here, we describe how to run rriScan on the testing datasets

- [download](#Downloading)
- [run rriScan](#Run)
- [Output](#Output)

# Downloading
Download the testing datasets to your local folder.

    ```bash
    # suppose your path to run rriScan is ./rriScan_test
    mkdir ./rriScan_test
    cd ./rriScan_test
    wget 'https://rnasysu.com/encori/software_test/rriScan.test.bam'
    wget 'https://rnasysu.com/encori/software_test/rriScan.test.bam.bai'
    ```

# Run
* It's very easy to run rriScan on testing dataset. Supposed `rriScan` has been added to you `PATH` environment and the `hg38.fa` and `hg38.fa.fai` have been generated.
    ```bash
    rriScan --fa hg38.fa --fai hg38.fa.fai --bam rriScan.test.bam --prefix rriScan --outdir ./
    ```

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

# Acknowledgements
Thanks a lot to everyone who contributed to the public codes and libraries (e.g. BamTools) used by rriScan.

# Contact
* Jian-Hua Yang <yangjh7@mail.sysu.edu.cn>, RNA Information Center, School of Life Sciences, Sun Yat-Sen University<BR>
* Keren Zhou <kzhou@coh.org>, Department of Systems Biology, Beckman Research Institute of City of Hope, Monrovia, CA, USA<BR>
