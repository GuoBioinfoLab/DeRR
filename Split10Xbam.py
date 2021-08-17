import pysam
import tqdm
import json
import argparse
import os
import logging

# Usage:
# python SplitVDJbam.py --bam /home/chensy/Dual/2.Result/scTCR-Seq/vdj_v1_cd8/outs/all_contig.bam --list /home/chensy/Dual/2.Result/scTCR-Seq/vdj_v1_cd8/outs/cell_barcodes.json --out ./

def CommandLineParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam",  required=True, help="BAM file")
    parser.add_argument("--list", required=True, help="Barcode list file")
    parser.add_argument("--out",  required=True, help="Output folder")
    parser.add_argument("--file", required=True, help="Generated Manifest file")
    return parser.parse_args()

if __name__ == "__main__":

    args = CommandLineParser()
    args = { arg:getattr(args, arg) for arg in vars(args) }


    LOG_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
    logging.basicConfig(level=logging.DEBUG, format=LOG_FORMAT)

    bam = args['bam']
    bc_liust = args['list']
    out_folder = os.path.abspath(args['out'])

    if not os.path.exists(out_folder):
        os.system(f"mkdir -p {out_folder}")

    logging.info("Loading Barcode list")

    with open(bc_liust) as handle:
        bc_list = set([ x.strip().split('-')[0] for x in handle.readlines() ])

    logging.info("START Writing")

    res = {}
    manifest = []

    idx = 0
    for rd in tqdm.tqdm(pysam.AlignmentFile(bam)):
        idx = idx + 1
        if rd.get_tag('CR') not in bc_list:
            continue
        bc, seq, qual = rd.get_tag("CR"), rd.seq, rd.qual
        manifest.append(bc)
        res.setdefault( bc, [] ).append(f"@{bc}_{idx}\n{seq}\n+\n{qual}\n")
        if len(res[bc]) > 1000:
            with open(f"{out_folder}/{bc}.fastq", 'a') as handle:
                handle.write("".join(res[bc]))
                res.pop(bc)
    for bc, vals in res.items():
        with open(f"{out_folder}/{bc}.fastq", 'a') as handle:
            handle.write("".join(vals))

