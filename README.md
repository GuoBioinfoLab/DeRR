# DeRR

<img src='./assets/DEERR_logo.png' align='right' height=350>

### What is **DeRR**

DeRR (Detection Dual T-cell receptor in single cell sequencing) is a method for:

- Detection TCR in single cell sequencing data
- Identification Dual-TCR in scRNA-Seq/scTCR-Seq

# Overview

- HomePage: [http://bioinfo.life.hust.edu.cn/DeRR]
- Github: [https://github.com/GuoBioinfoLab/DeRR]

# Installation

DeRR required following python packages:

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

Following tools are also required and should be able to access in $PATH or set their paths in `config.json`

* bwa
* samtools

We recommend using a virtual conda environment to install above packages and softwares:

```Shell
# 1. Create the environment and install the requirements
conda create -c conda-forge -c bioconda -n deer tqdm pandas biopython pysam networkx bwa samtools editdistance -y

# As sometimes conda might be very slow, users could use mamba instead of conda for faster installation
# install mamba by `conda install -n base -c conda-forge mamba`
# then run `mamba create -c conda-forge -c bioconda -n deer tqdm pandas biopython pysam networkx bwa samtools editdistance -y `

# 2. Download the IMGT reference and build the bwa index
# 2.1 Download the TRAV and TRBV gene sequences (fasta format, F+ORF+all P) from https://www.imgt.org/vquest/refseqh.html, save to refereces/TR-V.fa
# 2.2 Download the TRAJ and TRBJ gene sequences (fasta format, F+ORF+all P) from https://www.imgt.org/vquest/refseqh.html, save to refereces/TR-J.fa


# 3. Run DeRR
# (optional) Activate the envoriment `conda activate deer`
# Do some analysis
python DeRR.py --inf XXX --out XXX --threads number_of_threads
```

# Usage

> Before running the program, change the executable path of `bwa` and `samtools`  in `config.json`

Typical DeRR command for extraction Dual TCR will look like:

```Shell
python DeRR.py --inf /path/to/manifest.tsv --out /path/to/result.tsv --threads X
```

Users should list all the FASTQ files and Cell IDs (barcode) in the **manifest** file. The manifest file should contain at least 3 tab-separated columns like

| Data type | Cell-id    | Read1 File Name    | Read2 File Name  | Output Path (Option) |
|--------| -------- | -------- | --------------- | ------  |
| Paired  | Cell1  |  XXX   | XXX | XXX |
| Single  | Cell2  |  XXX | **None** | XXX |

A manifest file is like:

![](Manifest.png)

The **result.tsv** is like:

| v_call      | j_call     | junction_aa     | duplicate_count | locus | junction                                      | cell_id            | productive | sequence_id | sequence | rev_comp | d_call | sequence_alignment | germline_alignment | v_cigar | d_cigar | j_cigar |
|-------------|:-----------|:----------------|:----------------|:------|:----------------------------------------------|--------------------|------------|:------------|:---------|:---------|:-------|:-------------------|:-------------------|:--------|:--------|:--------|
| TRBV28*01   | TRBJ1-5*01 | CAITAGSANTEAFF  | 389             | TRB   | TGCATCGTCAGAGTCGCATCGGGTGGCGACTACAAGCTCAGCTTT | AAACCTGAGGCATTGG-1 | True       |             |          |          |        |                    |                    |         |         |         |
| TRAV26-1*02 | TRAJ20*01  | CIVRVASGGDYKLSF | 748             | TRA   | TGCATCGTCAGAGTCGCATCGGGTGGCGACTACAAGCTCAGCTTT | AAACCTGATTCATTGG-1 | True       |             |          |          |        |                    |                    |         |         |         |

To compatible with the AIRR format standard, there are some columns in the output of DeRR with **None** value like `d_call`, `d_cigar` etc. The columns with valuable information are list below:

| columns         | Description                                             |
|:----------------|:--------------------------------------------------------|
| v_call          | V allele gene name                                      |
| j_call          | J allele gene name                                      |
| junction_aa     | CDR3 amino acid sequence                                |
| duplicate_count | abundance of TCR in raw count                           |
| locus           | TRA/TRB                                                 |
| junction        | CDR3 nucleotide sequence                                |
| cell_id         | The name of input cell                                  |
| productive      | Is the TCR CDR3 region productive, which is always True |


## Tutorials

* [Handle 10X scTCR-Seq](./10XscTCR-Seq.ipynb)

## Cell demulpitexing


DeRR supports various types of single-cell data, but for sequencing types that do not directly provide a fastq file for each cell, users need to perform demultiplexing themselves to obtain fastq files corresponding to each cell, and then create  manifest.tsv file.

Taking 10X scRNA-Seq as an example, the Cell Barcode is stored in the CB tag of each record in the BAM file. Demultiplexing can be performed based on the CB tag of each record, we provide a script help demulpitexing the data:
```
python SplitVDJbam.py --bam all_contig.bam --list cell_barcodes.json --out /path/to/fastq_output --file /path/to/Manifest.tsv
```
where `all_contig.bam` and `cell_barcodes.json` is the output from cellranger, usually located in `{ProjectName}/outs`

We also have tested other types of scRNA-Seq data,

* Drop-Seq (data from PRJNA961723,  Cell barcode were inferred from the first 12 bases of read 1)
* CEL-Seq (PRJNA961718,  Cell barcode were inferred from the first 8 bases of read 1))
* MARS-Seq (PRJNA81671,  Cell barcode were inferred from the first 7 bases of read 2)


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

