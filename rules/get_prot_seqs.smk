from os.path import join
from pathlib import Path

input_dir = Path(config['input_dir'])
ASMS, = glob_wildcards(input_dir.joinpath("{asm}.faa"))
prefix = config['pref_mmseqs']

rule build_index:
    input: 
        expand(join(input_dir, "{asm}.faa"), asm=ASMS)
    output: 
        f"results/{prefix}.idx"
    run: 
        from Bio import SeqIO
        # Create a index with all the proteins sequences
        faas_idx = SeqIO.index_db(output[0], input, "fasta")


rule write_multi_fasta:
    input: rules.build_index.output
    output: "results/all_seqs.faa"
    run:
        from Bio import SeqIO
        faas_idx = SeqIO.index_db(input[0])
        records = (faas_idx[name] for name in faas_idx)
        with open(output[0], 'w') as fh:
            for rec in records:
                print(f'>{rec.description}\n{rec.seq}', file=fh)