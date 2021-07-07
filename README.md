# DeRR

<img src='DEERR_logo.png' align='right' height=350>

## What is **DERR**

DERR (Detection Dual T-cell receptor in single cell sequencing) is a toolkit for:

- Detection TCR for scRNA-Seq data
- Identification Dual-TCR  in scRNA-Seq
- Classify primary and minor chain for single CDR3 sequence

# Overview

- HomePage: [http://bioinfo.life.hust.edu.cn/DeRR]
- Github: [https://github.com/GuoBioinfoLab/DeRR]
- For detials please see our [Bioinformatics publication](!https://doi.org/10.1093/bioinformatics/btaa432)

# Installation

DeRR required follwing python packages:

* tqdm
* pandas 
* Biopython
* pysam
* networkx

Users could using `pip` or others package managers to install these packages like

```
pip install tqdm pandas biopython pysam networkx
```

Follwing tools are also required:

* bwa
* samtools
* fastp

We recommand using a vitrual conda envoriment to install above packages and softwares:

```Shell
# Create the envroiment and install the requirments
conda create -c conda-forge -c bioconda -n deer tqdm pandas biopython pysam networkx bwa samtools fastp -y

# As sometimes conda might be very slow, users could use mamba instead of conda for faster installation
conda install -n base -c conda-forge mamba #install mamba
mamba create -c conda-forge -c bioconda -n deer tqdm pandas biopython pysam networkx bwa samtools fastp -y #install requirments


# Activate the envoriment
conda activate deer
# Do some analysis
python DeRR.py --inf XXX --out XXX --threads 4
```


# Usage

Typical DERR command for extraction Dual TCR will look like:

```Shell
python DeRR.py --inf /path/to/manifest.tsv --out /path/to/result.tsv --threads X
```

Users should list all the FASTQ files and Cell IDs (barcode) in the manifest file. The manifest file should contain 3 tab-seprated columsn like

```
#For paired-end reads
Read1-file-name \t Read2-file-name \t Cell-id
#For single-end reads
Read1-file-name \t None \t Cell-id
```

For **10X V(D)J** sequencing data which don't provide FASTQ files for each cell, we provide a script help demulpitexing the data:

```
python SplitVDJbam.py --bam all_contig.bam --list cell_barcodes.json --out /path/to/fastq_output --file /path/to/Manifest.tsv
```
where `all_contig.bam` and `cell_barcodes.json` is the output from cellranger, usually located in `ProjectName/outs`

