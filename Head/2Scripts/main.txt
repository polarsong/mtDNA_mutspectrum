import os
import click

from tqdm import tqdm

from Bio import SeqIO

from src.feature_extraction.codon_usage import cds_codon_usage


@click.command
@click.argument("input",type=click.Path(exists=True), help="input multi genbank")
@click.argument("output",type=click.Path(), help="output directory")
def main(input,output):
    for seq in tqdm(SeqIO.parse(input, format='genbank')):
        seq_name = seq.name
        codons = cds_codon_usage(seq)
        codons.to_csv(os.path.join(output, '{}.csv'.format(seq_name)))


if __name__ == '__main__':
    main()