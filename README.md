# DeRR

<a href='blank'><img src='DEERR_logo.png' align='right' height='250'>

</a>



## What is **DERR

DERR (Detection Dual T-cell receptor in single cell sequencing) is a toolkit for:

- Detection TCR for scRNA-Seq data
- Identification Dual-TCR  in scRNA-Seq
- Classify primary and minor chain for single CDR3 sequence

# Overview

- HomePage: [http://bioinfo.life.hust.edu.cn/CATT](http://bioinfo.life.hust.edu.cn/CATT/Homepage.html)
- Github: [https://github.com/GuoBioinfoLab/CATT](https://github.com/GuoBioinfoLab/CATT)
- For detials please see our [Bioinformatics publication](!https://doi.org/10.1093/bioinformatics/btaa432)

# Installation

DeRR required follwing packages:

* tqdm
* pandas 
* Biopython
* pysam
* networkx

Users could using `pip` or others package managers to install these packages like

```
pip install tqdm pandas biopython pysam networkx
```



# Usage

Typical DERR command for extraction Dual TCR will look like:

```
python DERR.py --inf /path/to/manifest.tsv --out /path/to/result.tsv --threads X
```

Users should list all the FASTQ files and Cell IDS in the manifest file. The manifest file should contain 3 tab-seprated columsn like

```
#For paired-end reads
Read1-file-name \t Read2-file-name \t Cell-id
#For single-end reads
Read1-file-name \t None \t Cell-id
```

For data format like **10X V(D)J** that do not provide FASTQ files for each cell,  we provide a script to help demulpitexing the data:

```
python Split4bam.py --in --out 
```

