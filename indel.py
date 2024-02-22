## modules for generate indels

import pyhgvs as hgvs
import random
random.seed(404)
from pyhgvs.utils import read_transcripts
from pyfaidx import Fasta

# Read genome sequence using pyfaidx.
# Read RefSeq transcripts into a python dict.
def read_genome(genome_file, refgene_file):
    genome = Fasta(genome_file)
    with open(refgene_file) as infile:
        transcripts = read_transcripts(infile)
    return genome, transcripts

# Provide a callback for fetching a transcript by its name.
def get_transcript(name):
    return transcripts.get(name)

## get codon fasta
def index_containing_substring(the_list, substring):
    for i, s in enumerate(the_list):
        if substring in s:
              return i
    return -1

def get_transcript_fa(ref_fa, transcript):
    with open(ref_fa) as fi:
        line_list = fi.readlines()
        line_idx = index_containing_substring(line_list, transcript)
        print(line_list[line_idx])
        line_start ='a'
        coding_fa = ''
        while line_start != '>':
            coding_fa = coding_fa+ line_list[line_idx+1].strip()
            line_idx += 1
            line_start = line_list[line_idx+1][0]
    return coding_fa

def get_codons(coding_fa):
    with open(coding_fa) as fi:
        coding_seq = ''.join(''.join(fi.readlines()[1:]).split('\n'))
    return [coding_seq[i:i+3] for i in range(0,len(coding_seq), 3)]


# Parse the HGVS name into genomic coordinates and alleles.
## create HGVS names


def mk_cins(codon, ng):
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

def mk_cdel(codon, len_transcript, ng):
    '''
    yield c. del hgvs names
    '''
    i = 1
    while True:
        th = random.sample(range(1, len_transcript//3-3), 1)[0]
        start = th*3-1
        end = start +2
        del_seq = codon[th]
        del_hgvs = f'{ng}:c.{start}_{end}del{del_seq}'
        yield del_hgvs, i
        i += 1

## write out indels
def write_out_indels(out_vcf_file, num_ins, num_del, genome, get_transcript):
    with open(out_vcf_file, 'a') as out_vcf:
        out_vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for ins_names, num in mk_cins():
            if num > num_ins:
                break
            chrom, offset, ref, alt = hgvs.parse_hgvs_name(ins_names, genome, get_transcript=get_transcript, lazy=True)
            out_vcf.write(f'{chrom}\t{offset}\t{ins_names}\t{ref}\t{alt}\t.\tPASS\t.\n')
        for del_names, num in mk_cdel():
            if num > num_del:
                break
            chrom, offset, ref, alt = hgvs.parse_hgvs_name(del_names, genome, get_transcript=get_transcript, lazy=True)
            out_vcf.write(f'{chrom}\t{offset}\t{del_names}\t{ref}\t{alt}\t.\tPASS\t.\n')

## protein level

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

