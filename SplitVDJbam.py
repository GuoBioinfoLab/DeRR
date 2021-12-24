import pysam
import json
import argparse
import os
import logging

# Usage:
# python SplitVDJbam.py --bam /home/chensy/Dual/2.Result/scTCR-Seq/vdj_v1_cd8/outs/all_contig.bam --list /home/chensy/Dual/2.Result/scTCR-Seq/vdj_v1_cd8/outs/cell_barcodes.json --out ./

def CommandLineParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam",  required=True, help="BAM file")
    parser.add_argument("--list", required=True, help="Barcode list json")
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
        bc_list = set(json.load(handle))

    logging.info("Loading BAM file")

    rds =  [
        (rd.get_tag('CB'), rd.seq, rd.qual) for rd in pysam.AlignmentFile(bam) if rd.get_tag('CB') in bc_list
    ]

    logging.info("Sorting Cell Barcode")

    rds.sort(key= lambda x: x[0])

    logging.info("Writing to files")

    cur = 'START'
    manifest = []
    for idx, (bc, seq, qual) in enumerate(rds):

        if seq == 'None':
            continue

        if cur == 'START':
            cur = bc
            res = []
        

        if bc != cur:
            with open(f"{out_folder}/{cur}.fastq", 'w') as handle:
                handle.write("".join(res))
            manifest.append(f"{cur}\t{out_folder}/{cur}.fastq\tNone\n")
            res = []
            cur = bc

        res.append(f"@{bc}_{idx}\n{seq}\n+\n{qual}\n")
        
    else:
        #TODO(chens) there is a miss cell
        pass
    with open(args['file'], 'w') as handle:
        handle.write("".join(manifest))
