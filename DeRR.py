import argparse
from logging import debug
from Bio import SeqIO
import pandas as pd
import os
import time
import pysam
from multiprocessing import Pool
from tqdm.contrib.concurrent import process_map
import re
import networkx as nx
import string
import random
import editdistance
import sys
import pdb
import json
import numpy as np
import math

def selfLog(msg):
    print(time.ctime(time.time()) + "]     %s" % msg)

with open(os.path.realpath(sys.argv[0]).replace(os.path.split(sys.argv[0])[1], "config.json"), "r") as handle:
        global_config = json.load(handle)

global_config["TRV"] = os.path.abspath(sys.argv[0])[0:-7] + "reference/TR-V.fa"
global_config["TRJ"] = os.path.abspath(sys.argv[0])[0:-7] + "reference/TR-J.fa"

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

    def __init__(self, name, seq, vgene, jgene, cdr3, avlb):
        self.seq = seq
        self.name = name
        self.vgene = vgene
        self.jgene = jgene
        self.cdr3 = cdr3
        self.avlb = avlb

    def _repr__(self):
        return {
            'Seq': self.seq,
            'Name': self.name,
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

def hamming_distance(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

def real_score(rd, ref_seq):
    start, term = SegmentFromCigar(rd.cigartuples)
    lgt = term - start
    r_pos = rd.reference_start

    left_pad =  min(start, r_pos)
    right_pad = max(0, min(len(rd.seq) - term, len(ref_seq)-9 - (r_pos + lgt)))

    s1 = rd.seq[(start-left_pad):(term+right_pad)]
    s2 = ref_seq[(r_pos - left_pad):(r_pos + lgt+right_pad)] 

    return len(s1) - lgt - hamming_distance(s1, s2) * 3


def map2align(inp, ref, threads):

    bwa = global_config["bwa"]
    samtools = global_config["samtools"]

    prefix = os.path.realpath(sys.argv[0]).replace(os.path.split(sys.argv[0])[1], "")
    sam_file = prefix + "temporary/"+ ''.join(random.choices(string.ascii_uppercase + string.digits, k=10)) + '.tmp'
    bam_file = prefix + "temporary/" + ''.join(random.choices(string.ascii_uppercase + string.digits, k=10)) + '.tmp'
   
    os.system(f"{bwa} mem -t {threads} -SP -k 10 -A 1 -B 2 -L 0 -T 10 -v 0 {ref} {inp} 2>/dev/null > {sam_file}")
    os.system(f"{samtools} view -Sh -F 2308 {sam_file} 2>/dev/null > {bam_file}")
    os.system(f"rm -f {sam_file} &")
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

#Correct the false TCR within cell
def CellCorrect(group):

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
        group.iloc[idx, 4] = cnt + sum(group.loc[remove.index[1:], 'Counts'])


        if vgene == 'None':
            group.iloc[idx, 0] = most_common(remove['Vgene'], remove['Counts'])
        if jgene == 'None':
            group.iloc[idx, 1] = most_common(remove['Jgene'], remove['Counts'])

        group.drop(remove.index[1:], inplace=True)

    try:
        max_count = list(group['Counts'])[0]
        group = group[ (group.Vgene != 'None') & (group.Jgene != 'None') ]
        group = group[  group.Counts > max( max_count / 100, 2 )  ]
    except:
        pass

    return group

#Correct the false TCR based on cell population information
def PopCorrect(tab):

    #Background TCR level
    bg = []
    for bc, group in tab.groupby("CellId"):
        if group.shape[0] > 2:
            bg.extend( list(group['Counts'][2:]) )
    
    #asset no enough sample to filter
    try:
        tab = tab[ tab.Counts > np.percentile(bg, 90) ]
    except:
        return tab

    #Esimate the secondary TCR abundance
    thresold = np.percentile([ group.loc[group.index[1], 'Counts'] for bc, group in tab.groupby("CellId") if group.shape[0] > 1 ], 25)

    #Cell-within filter
    final = []
    for bc, group in tab.groupby("CellId"):
        if group.shape[0] > 2:
            #mean filter
            group = group[ group.Counts > math.sqrt(group.loc[group.index[0], 'Counts']*group.iloc[group.index[2], 'Counts']) ]
            if group.shape[0] > 2:
                group = group[ group.Counts > thresold ]
        final.append(group)
    
    return pd.concat(final)
    
def align(inp, threads, args):

    ##### QC
    r1, r2 = inp

    prefix = os.path.realpath(sys.argv[0]).replace(os.path.split(sys.argv[0])[1], "")
    os.system(f"mkdir -p {prefix}/temporary")
    output = prefix + "temporary/" + ''.join(random.choices(string.ascii_uppercase + string.digits, k=10)) + '.tmp'
    if r2 != "None":
        command = 'cat'
        if r1.split('.')[-1] == 'gz':
            command = 'zcat'
        os.system(f"{command} {r1}" + " | awk '{if(NR%4==1) $0=sprintf(\"@1_%d\",(1+i++)); print;}'" + f" > {output}")
        os.system(f"{command} {r2}" + " | awk '{if(NR%4==1) $0=sprintf(\"@2_%d\",(1+i++)); print;}'" + f" >> {output}")
    else:
        os.system(f"ln -s {os.path.realpath(r1)} {output}")

    with Pool(2) as pool:
        res = pool.starmap( map2align, [ (output, global_config["TRV"], threads // 2), (output, global_config["TRJ"], threads//2) ] )
    os.system(f'rm -f {output}')
    return res

def assignV(rd, refName2Seq):

    start, term = SegmentFromCigar(rd.cigartuples)
    refname = rd.reference_name
    tseq = rd.seq
    ref_seq = refName2Seq[ refname ]
    r_pos = rd.reference_start
    r_lgt = len(ref_seq)

    if (r_lgt - r_pos) - (len(tseq) - start  ) > 0:
        #discard
        return None
    elif real_score(rd, ref_seq) < -5:
        return None
    elif 'N' in tseq:
        return None
    else:
        return Myread(
            rd.qname,
            tseq[start:],
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

    if real_score(rd, ref_seq) < -6:
        return None
    elif 'N' in tseq:
        return None
    else:

         return Myread(
                rd.qname,
                tseq,
                "None",
                refname,
                "None",
                False
            )

def matchedNN(CDR3, nn):
    for idx, aa in enumerate(TranslateIntoAA(nn)):
        position = aa.find(CDR3)
        if position != -1:
            return nn[idx+position*3:idx+position*3+len(CDR3)*3]

#ll: list of pysam AligmentSegmetns
# def TongMing(ll, rev=1):

#     if len(ll) < 1:
#         return []

#     cur = ll[0].qname
#     res = []

#     for idx in range(0, len(ll)):

#         if ll[idx].qname != cur:
#             res.sort(key = lambda x: x.reference_start * rev)
#             yield res
#             cur = ll[idx].qname
#             res = [ll[idx]]
#         else:
#             res.append(ll[idx])
#     else:
#         yield res

def Extract_Motif(seq, cmotif, fmotif, coffset, foffset, innerC, innerF):

    res = 0
    Cx = [ m.end() - coffset for m in re.finditer(cmotif, seq) ]
    Fx = [ m.end() - foffset for m in re.finditer(fmotif, seq) ]

    if len(Cx)>1:
        Cx = [  pos for pos in Cx if pos in [ m.end() - 3 for m in re.finditer(innerC, seq) ] ]
    

    if len(Cx) ==0 and len(Fx) == 0:
        return ("None", 0)

    if (len(Cx) < 1 ) ^ (len(Fx) < 1):
        if len(Cx) < 1:
            Cx = [ m.end() - 3 for m in re.finditer(innerC, seq) ]
        else:
            Fx = [ m.end() + 2 for m in re.finditer(innerF, seq)]

  

    for (idx, xc) in enumerate(Cx):
        for xf in Fx:
            if (26 >=xf-xc-2 >= 7 ) and ( idx == len(Cx) -1 or not (26>=xf-Cx[idx+1]>=7)) and (not "*" in seq[xc:xf-2]):
                return (seq[xc:xf-2], 2)

    return ("None", 1)

def catt(inp, chain, threads):

    vbam, jbam = inp
    refName2Seq = {}

    # Read reference sequences
    for name in [ global_config["TRV"], global_config["TRJ"] ]:
        for seq in SeqIO.parse(name, 'fasta'):
            refName2Seq[ seq.id ] = str(seq.seq).upper()

    # Parse mapped reads from bam file and filter out unqulitifed reads
    vrs = []
    try:
        for rd in pysam.AlignmentFile(vbam, 'r'):
            if chain not in rd.reference_name:
                continue
            res = assignV(rd, refName2Seq)
            if res != None:
                vrs.append(res)
        vrs.sort(key = lambda x: x.name)
    except:
        print("No read in input file")


    jrs = []
    try:
        for rd in pysam.AlignmentFile(jbam, 'r'):
            if chain not in rd.reference_name:
                continue
            
            res = assignJ(rd, refName2Seq)
            if res != None:
                jrs.append(res)
        jrs.sort(key = lambda x: x.name)
    except:
        print("No read in input file")

    config = {
        'TRB':{
            "cmotif":"(LR|YF|YI|YL|YQ|YR)C(A|S|T|V|G|R|P|D)",
            "fmotif": "(F|T|Y|H)FG(A|D|E|N|P|Q|S)G",
            "coffset":2,
            "foffset":1,
            "innerC":"(CAS|CSA|CAW|CAT|CSV|CAI|CAR|CAV|CSG|CAN)",
            "innerF":"(QYF|YTF|AFF|QFF|LFF|QHF|LHF|LTF|IYF|KLF)",
        },
        'TRA': {
            "cmotif": "[EHILMSTV]{1}Y[FILY]{1}C[AGILV]{1}",
            "fmotif": "(LA|YI|FI|II|LY|LM|PT|TI|LV|ST|VT|LT|LI|LQ|MR|VI|FV|FQ|LF|LL|FE|FT|LS|LN|FY)F((ARG)|(G[A-Z]{1}G))",
            "coffset": 2,
            "foffset": 1,
            "innerC": "(CAF|CVF|CIF|CLF|CGF|CSF|CPF|CHF|CDF|CMF|CAV|CAA|CAY|CAE)",
            "innerF": "(LTF|LIF|LVF|FYF|LSF|YIF|LMF|FVF|IIF|TIF|AFF)"
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
                break

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
                    matchedNN(vrs[vdx].cdr3 if vrs[vdx].cdr3 != "None" else jrs[jdx].cdr3, vrs[vdx].seq if vrs[vdx].cdr3 != "None" else jrs[jdx].seq)
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

    # Add Edge
    for rd in vrs:
        if rd.cdr3 != "None":
            continue
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
        if rd.cdr3 != "None":
            continue
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

    #Find the forward flow
    flow_dict = nx.max_flow_min_cost(G, 'Source', 'Terminal')

    #Build the reverse network
    rG = G.reverse()
    for node, _ in j_nodes.items():
        rG['Terminal'][node]['capacity'] = 1
    for node, _ in v_nodes.items():
        rG[node]['Source']['capacity'] = 10241024
    #Find the reverse flow
    flow_dict_r = nx.max_flow_min_cost(rG, 'Terminal', 'Source')

    #Get reads that has flow in network
    left_v = [ x for (x, val) in flow_dict['Source'].items() if val > 0 ]
    left_j = [ x for (x, val) in flow_dict_r['Terminal'].items() if val > 0 ]


    for rd in vrs:
        if rd.name + '_V' not in left_v and rd.cdr3 != 'None' and rd.avlb:
            final_res.append((
                rd.vgene, 
                rd.cdr3, 
                'None',
                matchedNN(rd.cdr3, rd.seq)
            ))
    for rd in jrs:
        if rd.name + '_J' not in left_j and rd.cdr3 != 'None' and rd.avlb:
            final_res.append((
                'None', 
                rd.cdr3, 
                rd.jgene,
                matchedNN(rd.cdr3, rd.seq)
            ))

    repeat_list = []
    reduce_list = {}

    for rd in left_v:

        res = BuiltPath(FindPath(rd, 'Terminal', flow_dict)[0:-1], v_nodes, j_nodes)
        if res[0] < 2:

            for aa in TranslateIntoAA(res[1]):
                cdr3, code = Extract_Motif(aa, config[chain]['cmotif'], config[chain]['fmotif'], config[chain]['coffset'], config[chain]['foffset'], config[chain]['innerC'], config[chain]['innerF'])
                if code == 2:
                    #break
                    final_res.append((
                        res[2],
                        cdr3,
                        res[3],
                        matchedNN(cdr3, res[1])
                    ))
                    repeat_list.append(rd[0:-2])
                    break

    for rd in left_j:
        res = BuiltPath(FindPath( rd, 'Source', flow_dict_r )[0:-1][::-1], v_nodes, j_nodes)
        if res[0] < 2:
            for aa in TranslateIntoAA(res[1]):
                cdr3, code = Extract_Motif(aa, config[chain]['cmotif'], config[chain]['fmotif'], config[chain]['coffset'], config[chain]['foffset'], config[chain]['innerC'], config[chain]['innerF'])
                if code  == 2:
                    #break
                    final_res.append((
                        res[2],
                        cdr3,
                        res[3],
                        matchedNN(cdr3, res[1])
                    ))
                    if rd[0:-2] in repeat_list:
                        reduce_list[cdr3] = reduce_list.setdefault(cdr3, 0 ) + 1
                    break

    #Output results
    if len(final_res) > 0:

        tab = pd.DataFrame(final_res)
        tab.columns = [ 'Vgene', 'CDR3', 'Jgene', 'CDR3nn']
        # tab['Vgene'] = tab['Vgene'].apply(lambda x: x.split('*')[0])
        # tab['Jgene'] = tab['Jgene'].apply(lambda x: x.split('*')[0])
        tab['counts'] = 1
        unique_clone = set(tab['CDR3'])
        disjoint = { clone:clone for clone in unique_clone }
        for clone in unique_clone:
            for another in unique_clone - set([clone]):
                if clone in another:
                    disjoint[clone] = another
        tab = pd.DataFrame([ (most_common(group['Vgene'], group['counts']), 
                              most_common(group['Jgene'], group['counts']), 
                              disjoint[cdr3], 
                              most_common(group['CDR3nn'], group['counts']),
                              sum(group['counts']), 'TRB') for cdr3, group in tab.groupby('CDR3') ], columns = ['Vgene', 'Jgene', 'CDR3', 'CDR3nn', 'Counts', 'Chain'])
        for seq, val in reduce_list.items():
            tab.loc[ tab.CDR3 == seq, 'Counts' ] = tab.loc[ tab.CDR3 == seq, 'Counts' ] - val
        #tab = tab[ tab.Counts > 2 ]
        tab['Chain'] = chain
        return CellCorrect(tab.sort_values('Counts', ascending=False))
    else:
        return pd.DataFrame(columns = ['Vgene', 'Jgene', 'CDR3', 'Counts', 'Chain'])

def CommandLineParser():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--inf", help="Input file")
    group.add_argument("--r1", help="Read1 file")
    parser.add_argument("--r2", help="Read2 file", default="None")
    parser.add_argument("--out", help="Output folder", default="None")
    parser.add_argument("--align", type=int, default=4)
    parser.add_argument("--alleleseq", action="store_true")
    parser.add_argument("--population", action="store_true")
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args()

def Protocol(inp):

    if len(inp) < 6:
        r1, r2, sample_id, threads, args = inp
        output = None
    else:
        r1, r2, sample_id, threads, args, output = inp

    new_inp = align((r1, r2), threads, args)

    with Pool(2) as pool:
        res = pool.starmap( catt, [ (new_inp, 'TRA', threads), (new_inp, 'TRB', threads) ] )
    tmp = pd.concat(res)
    tmp['CellId'] = sample_id

    os.system(f"rm -f {new_inp[0]} {new_inp[1]} &")
    #print(new_inp)
    if output != None:
        tmp.to_csv(f"{output}", index=False, sep='\t')
    return tmp

def formatOutput(unform_res, args):


    unform_res.rename( columns= {
        "Vgene": 'v_call',
        "Jgene": 'j_call',
        "CDR3": "junction_aa",
        "Counts": 'duplicate_count',
        'Chain': "locus",
        'CDR3nn': "junction",
        'CellId': 'cell_id'
    }, inplace=True)

    unform_res['productive'] = 'True'

    for airr_col in ['sequence_id', 'sequence', 'rev_comp', 'd_call', 'sequence_alignment', 'germline_alignment', 'v_cigar', 'd_cigar', 'j_cigar']:
        unform_res[airr_col] = ""

    if args['alleleseq']:

        refName2Seq = {}
        for name in [ global_config["TRV"], global_config["TRJ"] ]:
            for seq in SeqIO.parse(name, 'fasta'):
                refName2Seq[ seq.id ] = str(seq.seq).upper()
        unform_res['v_sequence'] = unform_res['v_call'].apply(lambda x: refName2Seq[x])
        unform_res['j_sequence'] = unform_res['j_call'].apply(lambda x: refName2Seq[x])


if __name__ == "__main__":

    selfLog("Program start")

    # First run check

    for f in [  global_config['TRV'], global_config['TRJ'] ]:
        if not os.path.exists( f + '.bwt'):
            reads = [ x for x in SeqIO.parse(f, 'fasta') ]
            for read in reads:
                read.name = read.name.split("|")[1]
                read.id = read.id.split("|")[1]
                read.description = ""
            SeqIO.write(reads, f, 'fasta')
            os.system(f"{ global_config['bwa'] } index {f}")
    
    args = CommandLineParser()
    args = { arg:getattr(args, arg) for arg in vars(args) }

    max_thread  = min( args['threads'], args['align'])
    max_workers = max( args['threads'] // args['align'], 1 )

    refName2Seq = {}

    # Read reference sequences
    for name in [ global_config["TRV"], global_config["TRJ"] ]:
        for seq in SeqIO.parse(name, 'fasta'):
            refName2Seq[ seq.id ] = str(seq.seq).upper()

    #Single file input
    if args["r1"] != None:
        selfLog("Start detecting TCR")
        res = Protocol((args["r1"], args["r2"], "None", max_thread, args))
        if args["out"] != "None":
            formatOutput(res, args)
            res.to_csv(args["out"], index=False, sep='\t')
        selfLog("Detection end")
    else:
    #Manifest input 
        selfLog("Loading Manifest file")

        try:
            tab = pd.read_csv(f"{args['inf']}", sep='\t', index_col=0, header=None)
        except:
            selfLog("WARNING: Manifest file is empty")
            res = pd.DataFrame(columns = ["sequence_id", "sequence", "rev_comp", "productive", "v_call", "d_call", "j_call", "sequence_alignment," "germline_alignment", "junction",  "junction_aa", "v_cigar", "d_cigar", "j_cigar"])
            res.to_csv(args['out'], index=False, sep='\t')
            exit(0)

        selfLog("Start detecting TCR")
        if tab.shape[1] < 3:
            res = process_map(Protocol, [ (row[1], row[2], sample_id, max_thread, args) for sample_id, row in tab.iterrows() ], max_workers = max_workers, chunksize=4)
        else:
            res = process_map(Protocol, [ (row[1], row[2], sample_id, max_thread, args, row[3]) for sample_id, row in tab.iterrows() ], max_workers = max_workers, chunksize=4)
        selfLog("Detection end")

        if args["out"] != "None":
            tmp = pd.concat(res)
            if args['population']:
                tmp = pd.concat([ PopCorrect(tmp[tmp.Chain == 'TRA']), PopCorrect(tmp[tmp.Chain == 'TRB']) ])
            formatOutput(tmp, args)
            tmp.to_csv(args['out'], index=False, sep='\t')

    selfLog("Program end")