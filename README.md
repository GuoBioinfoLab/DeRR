# DeRR

<img src='DEERR_logo.png' align='right' height=350>

### What is **DERR**

DERR (Detection Dual T-cell receptor in single cell sequencing) is a toolkit for:

- Detection TCR in single cell sequencing data
- Identification Dual-TCR in scRNA-Seq/scTCR-Seq

# Overview

- HomePage: [http://bioinfo.life.hust.edu.cn/DeRR]
- Github: [https://github.com/GuoBioinfoLab/DeRR]

# Installation

DeRR required follwing python packages:

* tqdm
* pandas 
* Biopython
* pysam
* networkx
* editdistance

Users could using `pip` or others package managers to install these packages like

```
pip install tqdm pandas biopython pysam networkx editdistance
```

Follwing tools are also required and shoule be able to access in PATH:

* bwa
* samtools

We recommand using a vitrual conda envoriment to install above packages and softwares:

```Shell
# Create the envroiment and install the requirments
conda create -c conda-forge -c bioconda -n deer tqdm pandas biopython pysam networkx bwa samtools fastp editdistance -y

# As sometimes conda might be very slow, users could use mamba instead of conda for faster installation

#install mamba
conda install -n base -c conda-forge mamba
#install requirments
mamba create -c conda-forge -c bioconda -n deer tqdm pandas biopython pysam networkx bwa samtools fastp editdistance -y 


# Activate the envoriment
conda activate deer
# Do some analysis
python DeRR.py --inf XXX --out XXX --threads number_of_threads
```


# Usage

Before running the program, change the executable path of `bwa` and `samtools`  in `config.json`

Typical DERR command for extraction Dual TCR will look like:

```Shell
python DeRR.py --inf /path/to/manifest.tsv --out /path/to/result.tsv --threads X
```

Users should list all the FASTQ files and Cell IDs (barcode) in the **manifest** file. The manifest file should contain 3 tab-seprated columsn like

| Data type | Cell-id    | Read1 File Name    | Read2 File Name  | Output Path (Option) |
|--------| -------- | -------- | --------------- | ------  |
| Paired  | Cell1  |  XXX   | XXX | XXX |
| Single  | Cell2  |  XXX | **None** | XXX |

A manifest file is like:

![](Manifest.png)

The **result.tsv** is like:

| Vgene    | Jgene    | CDR3            | Counts | Chain | CellId             |
| -------- | -------- | --------------- | ------ | ----- | ------------------ |
| TRAV3    | TRAJ27   | CAHNTNAGKSTF    | 13     | TRA   | AAACCTGAGATCCTGT-1 |
| TRBV3-1  | TRBJ2-7  | CASSQGGALTYEQYF | 198    | TRB   | AAACCTGAGCGATAGC-1 |
| TRBV11-2 | TRBJ22-4 | CASSFDGLAKNIQYF | 68     | TRB   | AAACCTGAGGAGTCTG-1 |
| TRAV9-2  | TRAJ49   | CALFAGNQFYF     | 139    | TRA   | AAACCTGCATCTGGTA-1 |

For **10X V(D)J** sequencing data which don't provide FASTQ files for each cell, we provide a script help demulpitexing the data:

```
python SplitVDJbam.py --bam all_contig.bam --list cell_barcodes.json --out /path/to/fastq_output --file /path/to/Manifest.tsv
```
where `all_contig.bam` and `cell_barcodes.json` is the output from cellranger, usually located in `ProjectName/outs`

# Log

* 2022-09-15: Bug fixes

# Notice

The source code of DeRR include TCR V/J gene sequence from IMGT (https://www.imgt.org/vquest/refseqh.html), but for speeding up the program, all the sequences of V gene retain only the terminal 100bp.

# Term of Use

DeRR is maintained by An-Yuan Guo Bioinformatics Lab (Guo Lab). Guo Lab may, from time to time, update the content on https://github.com/GuoBioinfoLab/DeRR. Guo Lab makes no warranties or representations, express or implied, with respect to any of the Terms, including as to the present accuracy, completeness, timeliness, adequacy, or usefulness of any of the Terms. By using this website, you agree that Guo Lab will not be liable for any losses or damages arising from your use of or reliance on the Terms, or other websites or information to which this website may be linked.

DeRR is freely accessible for research use in an academic setting. You may view the Terms solely for your own personal reference or use for research in an academic setting. All academic research use of the Terms must credit DeRR as the source of the Terms and reference these Terms of Use; outside of scientific publication, you may not otherwise redistribute or share the Terms with any third party, in part or in whole, for any purpose, without the express permission of Guo Lab.

Unless you have signed a license agreement with Guo Lab, you may not use any part of the Terms for any other purpose, including:

* use or incorporation into a commercial product or towards performance of a commercial service;
* research use in a commercial setting;
* use for patient services; or
* generation of reports in a hospital or other patient care setting.

You may not copy, transfer, reproduce, modify or create derivative works of DeRR for any commercial purpose without the express permission of Guo Lab. If you seek to use Derr for such purposes, please request the license which best describes your anticipated use of DeRR below:

* Research use in commercial setting
* Use in a commercial product
* Use for patient services or reports in a hospital setting
* Please contact me at guoay@hust.edu.cn

