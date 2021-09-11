import numpy as np
import pandas as pd
from Bio import Seq
from Bio.Data import CodonTable
from Bio.Alphabet import  DNAAlphabet
from itertools import product
import re

dna_codons =[''.join(c) for c in product('ACGT', repeat = 3)]


def codon_usage(sequence: Seq, normalize = True):
    codon_list = list(CodonTable.standard_dna_table.forward_table.keys()) + CodonTable.standard_dna_table.stop_codons
    sequence = str(sequence)
    start_codon_positions = [re.search(startcodon, sequence) for startcodon in CodonTable.standard_dna_table.start_codons]
    start_codon_positions = list(filter(lambda x: x is not None, start_codon_positions))
    if len(start_codon_positions) == 0:
        raise ValueError('No start codons in seq')

    startpos = np.min(list(map(lambda x:x.start(),  start_codon_positions)))

    codons = []
    for p in  range(startpos, len(sequence), 3):
        codon = sequence[p:p+3]
        if codon not in dna_codons:
            continue

        codons.append(codon)
        if codon in CodonTable.standard_dna_table.stop_codons:
            break

    codons, counts = np.unique(codons, return_counts=True, axis=0)
    codons = [''.join(c) for c in codons]

    ret = pd.Series(index=codon_list, data=0)
    ret[codons] = counts
    if normalize:
        ret = ret / ret.sum()
    return ret

def cds_codon_usage(seq:Seq):
    cds_list = list(filter(lambda x: x.type == 'CDS' and x.strand == 1, seq.features))
    cds_seq_list = [cds.extract(seq.seq) for cds in cds_list]
    cds_stats = [codon_usage(c, normalize=False) for c in cds_seq_list]
    cds_stat_table = pd.concat(cds_stats, axis=1)
    cds_stat_table.columns = [cds.qualifiers['protein_id'][0] for cds in cds_list]
    return cds_stat_table


if __name__ == '__main__':
    from Bio import SeqIO

    seq = SeqIO.read('../../data/raw/genes.gbk', 'genbank')

    cds_list = list(filter(lambda x: x.type == 'CDS', seq.features))
    cds_seq_list = [cds.extract(seq) for cds in cds_list]
    cds_stats = [codon_usage(c, normalize=False) for c in cds_seq_list]
    joint_stat = sum(cds_stats)
    print(joint_stat.sum())
    #normed_stat = joint_stat / joint_stat.sum()
    cds_stat_table = pd.concat(cds_stats, axis=1)
    cds_stat_table.columns = [cds.qualifiers['protein_id'][0] for cds in cds_list]
    cds_stat_table.to_csv('../../data/interim/cds_codon_usage.csv')

