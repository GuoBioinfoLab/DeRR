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
    while True:
        cur = path[-1]
        if cur == term:
            return path
        for node, val in G[cur].items():
            if val > 0:
                path.append( node )
                break
        else:
            path.pop()
            if len(path) == 0:
                break

def BuiltPath(path, v_nodes, j_nodes, k=25):

    v = v_nodes[ path[0]  ][0]
    j = j_nodes[ path[-1] ][0]
    
    v_pos = v.find(path[1])
    j_pos = j.find(path[-2])
    segment = path[1] + "".join([ kmer[-1] for kmer in path[2:-1] ])
    
    total_len = (v_pos-j_pos+len(segment)-k) + len(j)
    fill_v = v + "*"*(total_len - len(v))
    fill_j = '*' * (v_pos-j_pos+len(segment)-k) + j
    
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