## modules for generate indels

import pyhgvs as hgvs
import random
random.seed(404)

# Provide a callback for fetching a transcript by its name.
def get_transcript(name):
    return transcripts.get(name)

def get_codons(coding_fa):
    with open(coding_fa) as fi:
        coding_seq = ''.join(''.join(fi.readlines()[1:]).split('\n'))
    return [coding_seq[i:i+3] for i in range(0,len(coding_seq), 3)]

# Parse the HGVS name into genomic coordinates and alleles.
## create HGVS names


def mk_cins(codon, ng, genome, transcripts):
    '''
    yield c. ins hgvs names
    '''
    i = 1
    while True:
        start = random.sample(range(0, 1089), 1)[0]
        end = start +1
        ins = random.sample(codon, 1)[0]
        ins_hgvs = f'{ng}:c.{start}_{end}ins{ins}'
        yield ins_hgvs, i
        i += 1

one2all ={'A': ('A', 'ALA', 'alanine'),
    'R': ('R', 'ARG', 'arginine'),
    'N': ('N', 'ASN', 'asparagine'),
    'D': ('D', 'ASP', 'aspartic acid'),
    'C': ('C', 'CYS', 'cysteine'),
    'Q': ('Q', 'GLN', 'glutamine'),
    'E': ('E', 'GLU', 'glutamic acid'),
    'G': ('G', 'GLY', 'glycine'),
    'H': ('H', 'HIS', 'histidine'),
    'I': ('I', 'ILE', 'isoleucine'),
    'L': ('L', 'LEU', 'leucine'),
    'K': ('K', 'LYS', 'lysine'),
    'M': ('M', 'MET', 'methionine'),
    'F': ('F', 'PHE', 'phenylalanine'),
    'P': ('P', 'PRO', 'proline'),
    'S': ('S', 'SER', 'serine'),
    'T': ('T', 'THR', 'threonine'),
    'W': ('W', 'TRP', 'tryptophan'),
    'Y': ('Y', 'TYR', 'tyrosine'),
    'V': ('V', 'VAL', 'valine'),
    'X': ('X', 'GLX', 'glutaminx'),
    'Z': ('Z', 'GLI', 'glycine'),
    'J': ('J', 'NLE', 'norleucine'),
    'U': ('U', 'CYC', 'cysteinc')
}

def mk_ins(fa, random_ins, ins):
    '''
    create p. ins hgvs names
    '''
    for i in random_ins:
        aa1 = one2all[fa[i]][1]
        aa2 = one2all[fa[i+1]][1]
        aai = random.choice(list(one2all.keys()))
        insi = f'p.{aa1}{i}_{aa2}{i+1}ins{aai}'
        ins.append(insi)
    return ins

def mk_del(fa, random_del, deletions):
    '''
    create p. del hgvs names
    '''
    for i in random_del:
        aa1 = one2all[fa[i]][1]
        deli = f'p.{aa1}{i}del'
        deletions.append(deli)
    return deletions