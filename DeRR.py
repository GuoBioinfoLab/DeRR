import argparse
from logging import debug
from Bio import SeqIO
import pandas as pd
import os
import time
import pysam
from tqdm.contrib.concurrent import process_map
import re
import networkx as nx
import string
import random
import editdistance
import sys

def selfLog(msg):
    print(time.ctime(time.time()) + "]     %s" % msg)

AAcode = {'TTT': 'F',
 'TTC': 'F',
 'TTA': 'L',
 'TTG': 'L',
 'TCT': 'S',
 'TCC': 'S',
 'TCA': 'S',
 'TCG': 'S',
 'TAT': 'Y',
 'TAC': 'Y',
 'TAA': '*',
 'TAG': '*',
 'TGT': 'C',
 'TGC': 'C',
 'TGA': '*',
 'TGG': 'W',
 'CTT': 'L',
 'CTC': 'L',
 'CTA': 'L',
 'CTG': 'L',
 'CCT': 'P',
 'CCC': 'P',
 'CCA': 'P',
 'CCG': 'P',
 'CAT': 'H',
 'CAC': 'H',
 'CAA': 'Q',
 'CAG': 'Q',
 'CGT': 'R',
 'CGC': 'R',
 'CGA': 'R',
 'CGG': 'R',
 'ATT': 'I',
 'ATC': 'I',
 'ATA': 'I',
 'ATG': 'M',
 'ACT': 'T',
 'ACC': 'T',
 'ACA': 'T',
 'ACG': 'T',
 'AAT': 'N',
 'AAC': 'N',
 'AAA': 'K',
 'AAG': 'K',
 'AGT': 'S',
 'AGC': 'S',
 'AGA': 'R',
 'AGG': 'R',
 'GTT': 'V',
 'GTC': 'V',
 'GTA': 'V',
 'GTG': 'V',
 'GCT': 'A',
 'GCC': 'A',
 'GCA': 'A',
 'GCG': 'A',
 'GAT': 'D',
 'GAC': 'D',
 'GAA': 'E',
 'GAG': 'E',
 'GGT': 'G',
 'GGC': 'G',
 'GGA': 'G',
 'GGG': 'G'}

def BreakSeqIntoAAcodes(seq, frame, n):

    return ''.join([AAcode[seq[i:i+3]] for i in range(frame, n, 3)])

def TranslateIntoAA(seq):

    n = len(seq)
    return [ BreakSeqIntoAAcodes(seq, ff, n-(n-ff)%3) for ff in [0, 1, 2] ]

#There must be some more elegant and faster way to do this
#Tip1: according to MCMF value to see wheter the network is liantong 
def FindPath(cur, term, G):
    path = [cur]
    marker = set()
    while True:
        cur = path[-1]
        if cur == term:
            return path
        for node, val in G[cur].items():
            if val > 0 and node not in path and (cur, node) not in marker:
                marker.add( (cur, node) )
                path.append( node )
                break
        else:
            path.pop()
            if len(path) == 0:
                break

def BuiltPath(path, v_nodes, j_nodes, k=25):

# Some code for these is no path
    if path == None or len(path) < 2:
        return (3, "", "", "")

    v = v_nodes[ path[0]  ][0]
    j = j_nodes[ path[-1] ][0]

    v_pos = v.find(path[1])
    j_pos = j.find(path[-2])
    segment = path[1] + "".join([ kmer[-1] for kmer in path[2:-1] ])

    total_len = (v_pos-j_pos+len(segment)-k) + len(j)
    fill_v = v[0:v_pos] + segment + '*'*(total_len - len(segment) - v_pos)
    fill_j = '*'*(total_len - len(j)) + j

    res = 0
    merged_seq = []
    for idx in range(len(fill_v)):
        merged_seq.append( fill_v[idx] if fill_v[idx] != '*' else fill_j[idx] )
        if fill_v[idx] == '*' or fill_j[idx] == '*':
            continue
        else:
            if fill_v[idx] != fill_j[idx]:
                res = res + 1
                if res > 1:
                    break


    return (res, "".join(merged_seq), v_nodes[ path[0] ][1], j_nodes[ path[-1] ][1])


class Myread(object):

    def __init__(self, name, seq, qual, vgene, jgene, cdr3, avlb):
        self.seq = seq
        self.name = name
        self.qual = qual
        self.vgene = vgene
        self.jgene = jgene
        self.cdr3 = cdr3
        self.avlb = avlb

    def _repr__(self):
        return {
            'Seq': self.seq,
            'Name': self.name,
            'Qual': self.qual,
            'Gene': " ".join([ self.vgene, self.jgene ]),
            'CDR3': self.cdr3,
            'Aviable': self.avlb
        }

    def __len__(self):
        return len(self.seq)

def SegmentFromCigar( cigartuples  ):
    start = 0
    end = 0
    pos = 0
    flag = True
    for it in cigartuples:
        if flag and it[0] == 0:
            start = pos
            flag = False
        pos = pos + it[1]
        if it[0] == 0:
            end = pos
    return start,end

def real_score(rd, ref_seq):
    start, _ = SegmentFromCigar(rd.cigartuples)
    r_pos = rd.reference_start
    score = 0
    for index in range(len(rd.seq)):
        offset = index - start
        if 0 <= r_pos + offset < len(ref_seq) - 9 :
            if rd.seq[index] == ref_seq[r_pos + offset]:
                score = score + 1
            else:
                score = score - 2
    return score

def map2align(inp, ref, threads):

    #TODO(chensy) change to config.json
    bwa = "/workspace/chensy/dual/0.Script/deer/bwa-mem2-2.2.1_x64-linux/bwa-mem2"
    samtools = "samtools"

    prefix = os.path.realpath(sys.argv[0]).replace("DeRR.py", "")
    sam_file = prefix + "temporary/"+ ''.join(random.choices(string.ascii_uppercase + string.digits, k=10)) + '.tmp'
    bam_file = prefix + "temporary/" + ''.join(random.choices(string.ascii_uppercase + string.digits, k=10)) + '.tmp'
    #os.system(f"{bwa} mem -t {threads} -r 2.3 -k 10 -A 1 -B 2 -L 0 -T 17 -v 0 {ref} {inp} 2>/dev/null > {sam_file}")
    os.system(f"{bwa} mem -t {threads} -r 2.3 -k 10 -A 1 -B 2 -L 0 -T 17 -v 0 {ref} {inp} 2>/dev/null > {sam_file}")
    os.system(f"{samtools} view -Sh -F 2308 {sam_file} 2>/dev/null > {bam_file}")
    os.system(f"rm -f {sam_file}")
    return bam_file

def most_common(vals, cnts):
    total = {}
    for (val, cnt) in zip(vals, cnts):
        total[val] = total.setdefault(val, 0) + cnt
    total = list(total.items())
    total.sort(key=lambda x:x[1], reverse=True)
    for item, val in total:
        if item!= 'None':
            return item
    else:
        return 'None'

def Correct(group):

    for idx in range(0, 100):

        if idx >= group.shape[0]:
            continue

        seq   = list(group['CDR3'])[idx]
        cnt   = list(group['Counts'])[idx]
        vgene = list(group['Vgene'])[idx]
        jgene = list(group['Jgene'])[idx]

        remove = group[ (group['CDR3'].apply(lambda x: editdistance.eval(x, seq) <= len(seq) // 6 )) ]

        if remove.shape[0] < 2:
            continue

        #Counts
        group.iloc[idx, 3] = cnt + sum(group.loc[remove.index[1:], 'Counts'])
        if vgene == 'None':
            group.iloc[idx, 0] = most_common(remove['Vgene'], remove['Counts'])
        if jgene == 'None':
            group.iloc[idx, 0] = most_common(remove['Jgene'], remove['Counts'])

        group.drop(remove.index[1:], inplace=True)

    return group

def align(inp, threads=1):

    ##### QC
    fastp = 'fastp'
    r1, r2 = inp

    prefix = os.path.realpath(sys.argv[0]).replace("DeRR.py", "")
    os.system(f"mkdir -p {prefix}/temporary")
    output = prefix + "temporary/" + ''.join(random.choices(string.ascii_uppercase + string.digits, k=10)) + '.tmp'
    if r2 != "None":
        os.system(f"{fastp} -i {r1} -I {r2} -m --merged_out {output} --include_unmerged --detect_adapter_for_pe -q 20 -e 25 -L 30 -c -g -w {threads} -h /dev/null -j /dev/null --overlap_len_require 20 --overlap_diff_limit 2  >/dev/null 2>&1")
    else:
        os.system(f"{fastp} -i {r1} -o {output} -q 20 -e 25 -L 30 -c -g -w {threads} -h /dev/null -j /dev/null  >/dev/null 2>&1")


    res = (
        map2align(output, f"{prefix}reference/MEM2/AIRR-V-DNA.fa", threads),
        map2align(output, f"{prefix}reference/MEM2/AIRR-J-DNA.fa", threads)
    )
    os.system(f'rm -f {output}')
    return res


def assignV(rd, refName2Seq):

    start, term = SegmentFromCigar(rd.cigartuples)
    refname = rd.reference_name
    tseq = rd.seq
    ref_seq = refName2Seq[ refname ]
    r_pos = rd.reference_start
    r_lgt = len(ref_seq)


    if (r_lgt - r_pos) - (len(tseq) - start  ) > -10:
        #discard
        return None
    elif real_score(rd, ref_seq) < 0:
        return None
    elif 'N' in tseq:
        return None
    else:
        return Myread(
            rd.qname,
            tseq[start:],
            rd.qual[start:],
            refname,
            "None",
            "None",
            False
        )

def assignJ(rd, refName2Seq):

    start, term = SegmentFromCigar(rd.cigartuples)
    refname = rd.reference_name
    tseq = rd.seq
    ref_seq = refName2Seq[ refname ]
    r_pos = rd.reference_start
    r_lgt = len(ref_seq)

    if start - 1 < 10:
        return None
    elif real_score(rd, ref_seq) < 0:
        return None
    elif 'N' in tseq:
        return None
    else:
#         return Myread(
#             rd.qname,
#             tseq[0:start] + ref_seq[r_pos:],
#             rd.qual[0:start] + 'G'*(r_lgt-r_pos),
#             "None",
#             refname,
#             "None",
#             False
#         )

         return Myread(
                rd.qname,
                tseq,
                rd.qual,
                "None",
                refname,
                "None",
                False
            )



#ll: list of pysam AligmentSegmetns
def TongMing(ll, rev=1):

    if len(ll) < 1:
        return []

    cur = ll[0].qname
    res = []

    for idx in range(0, len(ll)):

        if ll[idx].qname != cur:
            res.sort(key = lambda x: x.reference_start * rev)
            yield res
            cur = ll[idx].qname
            res = [ll[idx]]
        else:
            res.append(ll[idx])
    else:
        yield res

#TODO(chensy) 2 round to found potential reads
#try to avoid error from sequences that have only innerC
def Extract_Motif(seq, cmotif, fmotif, coffset, foffset, innerC, innerF):

    res = 0
    Cx = [ m.end() - coffset for m in re.finditer(cmotif, seq) ]
    Fx = [ m.end() - foffset for m in re.finditer(fmotif, seq) ]

    if len(Cx) ==0 and len(Fx) == 0:
        return ("None", 0)

    # if (len(Cx) < 0 ) ^ (len(Fx) < 0):
    #     if len(Cx) < 0:
    #         Cx = [ m.end() -2 for m in re.finditer(innerC, seq) ]
    #     else:
    #         Fx = [ m for m in re.finditer(innerF, seq)]

    for (idx, xc) in enumerate(Cx):
        for xf in Fx:
            if (22 >=xf -xc >= 6 ) and ( idx == len(Cx) -1  or not (32>=xf-Cx[idx+1]>=7)) and not "*" in seq[xc:xf]:
                return (seq[xc:xf-2], 2)

    return ("None", 1)

def Extract_CDR3(seq, cmotif, fmotif, coffset, foffset, innerC, innerF):

    res = 0
    Cx = [ m.end() - coffset for m in re.finditer(cmotif, seq) ]
    Fx = [ m.end() - foffset for m in re.finditer(fmotif, seq) ]

    if len(Cx) ==0 and len(Fx) == 0:
        return ("None", 0)

    if (len(Cx) < 0 ) ^ (len(Fx) < 0):
        if len(Cx) < 0:
            Cx = [ m.end() -2 for m in re.finditer(innerC, seq) ]
        else:
            Fx = [ m for m in re.finditer(innerF, seq)]

    for (idx, xc) in enumerate(Cx):
        for xf in Fx:
            if (22 >=xf -xc >= 6 ) and ( idx == len(Cx) -1  or not (32>=xf-Cx[idx+1]>=7)) and not "*" in seq[xc:xf]:
                return (seq[xc:xf-2], 2)

    return ("None", 1)


def catt(inp, chain, threads):

    vbam, jbam = inp

    refName2Seq = {}
    prefix = os.path.realpath(sys.argv[0]).replace("DeRR.py", "")
    for name in [ f"{prefix}reference/AIRR-V-DNA.fa", f"{prefix}reference/AIRR-J-DNA.fa"]:
        for seq in SeqIO.parse(name, 'fasta'):
            refName2Seq[ seq.id ] = str(seq.seq).upper()


    try:
        vrs = [ rd for rd in pysam.AlignmentFile(vbam, 'r') if chain in rd.reference_name ]
    except:
         vrs = []
    vrs.sort(key = lambda x: x.qname)
    vrs = [ x[0] for x in TongMing(vrs, -1) ]
    vrs = list(filter(lambda x: x!= None, map(lambda x: assignV(x, refName2Seq), vrs)))


    try:
        jrs = [ rd for rd in pysam.AlignmentFile(jbam, 'r') if chain in rd.reference_name ]
    except:
        jrs = []
    jrs.sort(key = lambda x: x.qname)
    jrs = [ x[0] for x in TongMing(jrs) ]
    #remove reads taht have no contribution to the result
    jrs = list(filter(lambda x: x!= None, map(lambda x: assignJ(x, refName2Seq), jrs)))



    config = {
        'TRB':{
            "cmotif":"(LR|YF|YI|YL|YQ|YR)C(A|S|T|V|G|R|P|D)",
            "fmotif": "(F|T|Y|H)FG(A|D|E|N|P|Q|S)G",
            "coffset":2,
            "foffset":1,
            "innerC":"place_holder",
            "innerF":"place_holder",
        },
        'TRA': {
            "cmotif": "[EHILMSTV]{1}Y[FILY]{1}C[AGILV]{1}",
            "fmotif": "(LA|YI|FI|II|LY|LM|PT|TI|LV|ST|VT|LT|LI|LQ|MR|VI|FV|FQ|LF|LL|FE|FT|LS|LN|FY)F((ARG)|(G[A-Z]{1}G))",
            "coffset": 2,
            "foffset": 1,
            "innerC": "place_holder",
            "innerF": "place_holder"
        }

    }


    #### 2. First retrieve
    ####    Find reads that have complete fregments

    for rd in vrs:
        for aa in TranslateIntoAA(rd.seq):
            res = Extract_Motif(aa, config[chain]['cmotif'], config[chain]['fmotif'], config[chain]['coffset'], config[chain]['foffset'], config[chain]['innerC'], config[chain]['innerF'])
            if res[1] > 0:
                rd.avlb = True
            if res[1] == 2:
                rd.cdr3 = res[0]

    for rd in jrs:
        for aa in TranslateIntoAA(rd.seq):
            res = Extract_Motif(aa, config[chain]['cmotif'], config[chain]['fmotif'], config[chain]['coffset'], config[chain]['foffset'], config[chain]['innerC'], config[chain]['innerF'])
            if res[1] > 0:
                rd.avlb = True
            if res[1] == 2:
                rd.cdr3 = res[0]
                break




    ##### Build kmer table
    cnt = {}
    k = 25
    kmer_offset = 5

    for rd in filter(lambda x: x.avlb and len(x.seq) >= 35, vrs):
        seq = rd.seq
        l = len(seq)
        for i in range( l-k-kmer_offset, l-k+1 ):
            kmer = seq[i:(i+k)]
            cnt[kmer] = cnt.setdefault(kmer, 0) + 1

    for rd in filter(lambda x: x.avlb and len(x.seq) >= 35, jrs):
        seq = rd.seq
        l = len(seq)
        for i in range( 0, min(kmer_offset + k, l)):
            kmer = seq[i:(i+k)]
            cnt[kmer] = cnt.setdefault(kmer, 0) + 1


##### Remove reads that span over junction
    vdx, jdx = 0, 0
    final_res = []
    while(vdx < len(vrs) and jdx < len(jrs)):
        if vrs[vdx].name == jrs[jdx].name:

            if vrs[vdx].cdr3 != "None" or jrs[jdx].cdr3 != "None":
                final_res.append((
                    vrs[vdx].vgene,
                    vrs[vdx].cdr3 if vrs[vdx].cdr3 != "None" else jrs[jdx].cdr3,
                    jrs[jdx].jgene,
                ))
                vrs[vdx].avlb = False
                jrs[jdx].avlb = False

            vdx = vdx + 1
            jdx = jdx + 1
        elif vrs[vdx].name < jrs[jdx].name:
            while(vdx < len(vrs) and vrs[vdx].name < jrs[jdx].name ):
                vdx = vdx + 1
        else:
            while(jdx < len(jrs) and vrs[vdx].name > jrs[jdx].name ):
                jdx = jdx + 1
    #remove
    vrs = list(filter(lambda x: x.avlb and len(x.seq) >= 35, vrs))
    jrs = list(filter(lambda x: x.avlb and len(x.seq) >= 35, jrs))


    ####### Build the network, and run MFNC
    G = nx.DiGraph()
    G.add_node('Source', weight = 1)
    G.add_node("Terminal", weight =1 )

    v_nodes = {}
    j_nodes = {}

    for rd in vrs:
        seq = rd.seq
        l = len(seq)
        rd_node = rd.name + '_V'
        v_nodes[rd_node] = (rd.seq, rd.vgene)
        nodes = [ seq[i:(i+k)] for i in range(l-k-kmer_offset, l-k+1) ]
        G.add_node(rd_node)
        G.add_edge('Source', rd_node, capacity=1, weight=1)
        for node in nodes:
            G.add_node(node)
            G.add_edge(rd_node, node, capacity=1, weight = -cnt[node])
        nodes = iter(nodes)
        for x, y in zip(nodes, nodes):
            G.add_edge(x, y, capacity = 102410241024, weight = -cnt[y])

    for rd in jrs:
        seq = rd.seq
        l = len(seq)
        rd_node = rd.name + '_J'
        j_nodes[rd_node] = (rd.seq, rd.jgene)
        nodes = [ seq[i:(i+k)] for i in range( 0, min(kmer_offset, l-k+1)) ]
        G.add_node(rd_node)
        G.add_edge(rd_node, 'Terminal',  capacity = 10241024, weight=1)
        for node in nodes:
            G.add_node(node)
            G.add_edge(node, rd_node, capacity=1, weight= -cnt[node])
        nodes = iter(nodes)
        for x, y in zip(nodes, nodes):
            G.add_edge(x, y, capacity = 102410241024, weight = -cnt[y])


    flow_dict = nx.max_flow_min_cost(G, 'Source', 'Terminal')
    rG = G.reverse()
    for node, _ in j_nodes.items():
        rG['Terminal'][node]['capacity'] = 1
    for node, _ in v_nodes.items():
        rG[node]['Source']['capacity'] = 10241024
    flow_dict_r = nx.max_flow_min_cost(rG, 'Terminal', 'Source')

    left_v = [ x for (x, val) in flow_dict['Source'].items() if val > 0 ]
    left_j = [ x for (x, val) in flow_dict_r['Terminal'].items() if val > 0 ]

    #TODO(chensy) Different of Extract_Motif and Extract_CDR3
    for rd in vrs:
        if rd.name + '_V' not in left_v and rd.cdr3 != 'None':
            final_res.append((rd.vgene, rd.cdr3, 'None'))
    for rd in jrs:
        if rd.name + '_J' not in left_j and rd.cdr3 != 'None':
            final_res.append(('None', rd.cdr3, rd.jgene))


    repeat_list = []
    reduce_list = {}

    for rd in left_v:

        res = BuiltPath(FindPath(rd, 'Terminal', flow_dict)[0:-1], v_nodes, j_nodes)
        if res[0] < 2:

            for aa in TranslateIntoAA(res[1]):
                cdr3, code = Extract_Motif(aa, config[chain]['cmotif'], config[chain]['fmotif'], config[chain]['coffset'], config[chain]['foffset'], config[chain]['innerC'], config[chain]['innerF'])
                if code == 2:
                    final_res.append((
                        res[2],
                        cdr3,
                        res[3]
                    ))
                    repeat_list.append(rd[0:-2])
                    break

    for rd in left_j:
        res = BuiltPath(FindPath( rd, 'Source', flow_dict_r )[0:-1][::-1], v_nodes, j_nodes)
        if res[0] < 2:
            for aa in TranslateIntoAA(res[1]):
                cdr3, code = Extract_Motif(aa, config[chain]['cmotif'], config[chain]['fmotif'], config[chain]['coffset'], config[chain]['foffset'], config[chain]['innerC'], config[chain]['innerF'])
                if code  == 2:
                    final_res.append((
                        res[2],
                        cdr3,
                        res[3]
                    ))
                    if rd[0:-2] in repeat_list:
                        reduce_list[cdr3] = reduce_list.setdefault(cdr3, 0 ) + 1
                    break


    if len(final_res) > 0:

        tab = pd.DataFrame(final_res)
        tab.columns = [ 'Vgene', 'CDR3', 'Jgene']
        tab['Vgene'] = tab['Vgene'].apply(lambda x: x.split('*')[0])
        tab['Jgene'] = tab['Jgene'].apply(lambda x: x.split('*')[0])
        tab['counts'] = 1

        tab = pd.DataFrame([ (most_common(group['Vgene'], group['counts']), most_common(group['Jgene'], group['counts']), cdr3, sum(group['counts']), 'TRB') for cdr3, group in tab.groupby('CDR3') ], columns = ['Vgene', 'Jgene', 'CDR3', 'Counts', 'Chain'])
        for seq, val in reduce_list.items():
            tab.loc[ tab.CDR3 == seq, 'Counts' ] = tab.loc[ tab.CDR3 == seq, 'Counts' ] - val
        tab = tab[ tab.Counts > 2 ]
        tab['Chain'] = chain

        return Correct(tab.sort_values('Counts', ascending=False))
    else:
        return pd.DataFrame(columns = ['Vgene', 'Jgene', 'CDR3', 'Counts', 'Chain'])

def CommandLineParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inf",  required=True, help="Input file")
    parser.add_argument("--out",  required=True, help="Output folder")
    parser.add_argument("--align", type=int, default=4)
    parser.add_argument("--threads", type=int, default=2)
    return parser.parse_args()

def Protocol(inp):

    if len(inp) < 5:
        r1, r2, sample_id, threads = inp
        output = None
    else:
        r1, r2, sample_id, threads, output = inp

    new_inp = align((r1, r2),      threads)
    alpha   = catt(new_inp, 'TRA', threads)
    beta    = catt(new_inp, 'TRB', threads)

    tmp = pd.concat([alpha, beta])
    tmp['CellId'] = sample_id

    os.system(f"rm -f {new_inp[0]} {new_inp[1]}")
    if output != None:
        tmp.to_csv(f"{output}", index=False, sep='\t')
    return tmp


if __name__ == "__main__":

    selfLog("Program start")

    args = CommandLineParser()
    args = { arg:getattr(args, arg) for arg in vars(args) }

    selfLog("Loading Manifest file")

    try:
        tab = pd.read_csv(f"{args['inf']}", sep='\t', index_col=0, header=None)
    except:
        selfLog("WARNING: Manifest file is empty")
        res = pd.DataFrame(columns = ['Vgene', 'Jgene', 'CDR3', 'Counts', 'Chain'])
        res.to_csv(args['out'], index=False, sep='\t')
        exit(0)


    max_thread  = min( args['threads'], args['align'])
    max_workers = max( args['threads'] // args['align'], 1 )

    selfLog("Start detecting TCR")
    if tab.shape[1] < 3:
        res = process_map(Protocol, [ (row[1], row[2], sample_id, max_thread) for sample_id, row in tab.iterrows() ], max_workers = max_workers, chunksize=2)
    else:
        res = process_map(Protocol, [ (row[1], row[2], sample_id, max_thread, row[3]) for sample_id, row in tab.iterrows() ], max_workers = max_workers, chunksize=2)
    selfLog("Detection end")


    pd.concat(res).to_csv(args['out'], index=False, sep='\t')
    selfLog("Program end")

