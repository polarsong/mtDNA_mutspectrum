import pandas as pd
from src.feature_extraction.codon_usage import cds_codon_usage

from Bio import SeqIO
from tqdm import tqdm
import io


def get_seqs(fn):
    src_fn = open(fn)
    file_text = src_fn.read()
    start, end = 0, file_text.find('//\n')

    recs = []
    with tqdm(total=len(file_text)) as t:
        while end != -1:
            try:
                rec = file_text[start: end + 3]
                start = end
                end = file_text.find('//\n', end + 1)

                s = io.StringIO()
                s.write(rec)
                s.seek(0)
                recs.append(SeqIO.read(s, format='genbank'))
            except:
                pass
            finally:
                t.update(end)
    return recs

# ftp://ftp.ncbi.nih.gov/refseq/release/viral/
if __name__ == '__main__':
    # with open('../../data/raw/split_viral_bgff/1/NC_030454.gbk') as f:
    #     sequence = SeqIO.read(f, 'genbank')
    #     cds_codon_usage(sequence)
    #
    # exit()
    src_file = '../../data/raw/viral.1.genomic.gbff'
    seqs = get_seqs(src_file)
    # for seq in tqdm(seqs):
    #     with open('../../data/raw/split_viral_bgff/2/{}.gbk'.format(seq.name), 'w') as f:
    #         SeqIO.write(seq, f, 'genbank')
    #

    seq_list = pd.read_csv('../../data/raw/sequences.csv')
    for seq in tqdm(SeqIO.parse(src_file, format='genbank')):
        seq_name = seq.name
        try:
            if seq_name in seq_list.Accession.values:
                codons = cds_codon_usage(seq)
                codons.to_csv('../../data/interim/1/{}.csv'.format(seq_name))
        except Exception as e:
            print(seq_name)
            print(e)
            exit()

    src_file = '../../data/raw/viral.2.genomic.gbff'
    seqs = get_seqs(src_file)
    for seq in tqdm(seqs):
        seq_name = seq.name
        if seq_name in seq_list.Accession.values:
            codons = cds_codon_usage(seq)
            codons.to_csv('../../data/interim/2/{}.csv'.format(seq_name))
