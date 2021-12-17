import sys

sys.stderr = open(snakemake.log[0], "w")

from Bio import SeqIO

with open(snakemake.input[0], 'r') as fh:
    genes = [gene.strip() for gene in fh]

with open(snakemake.output[0], 'a') as fh:
    for rec in SeqIO.parse(snakemake.input[1], format= 'fasta'):
        if rec.name in genes:
            print(f'>{rec.name}\n{rec.seq}', file=fh)
