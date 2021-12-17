with open(snakemake.input[0], 'r') as fh:
    lines = [line.strip() for line in fh.readlines()]

kos = [ko.split(':')[1] for ko in lines for ko in ko.split('\t')[1].split(',')]

with open(snakemake.output[0], 'a') as fh:
    for ko in kos:
        print(ko, file=fh)
