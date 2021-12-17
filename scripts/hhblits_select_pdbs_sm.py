import sys 

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path
from collections import defaultdict
import pickle
import re


def parse_hhblits_table(file):
    
    with open(file, 'r') as fh:
        vals = []
        switch = False
        counter = 0
        for line in fh:
            if re.search(r'^\s+No Hit', line):
                switch = True
            elif switch:
#                 print(line)
                counter += 1
                vals.append(re.search(r'[0-9]\s+(\w+).+?\s+([0-9]+\.?[0-9]+)', line.strip()).groups())
#                 print(re.split('\s+', line))
            if counter == 3:
                break
    return vals


def select_pdbs(values, f, prob_tresh):
    '''
    values: A list of tuples coming from parse_hhblits_table()
    f: file to be parsed
    prob_tresh: probability value in order to select a pdb (int)
    '''
    high_probs = defaultdict(list)
    for tup in values:
        prob = float(tup[1])
        pdb_id = tup[0].split('_')[0]
        # print(Path(f).stem, prob)
        if prob >= prob_tresh:
            # print(pdb_id, prob)
            gene = '_'.join(Path(f).stem.split('_')[:2])
            high_probs[gene].append((pdb_id, prob))
            break
    
    return high_probs

prob = snakemake.params.prob
files = snakemake.input
print(files, prob)

selected_high = []
for f in files:
    values = parse_hhblits_table(f)
    high_probs = select_pdbs(values, f, prob)
    if high_probs:
        selected_high.append(high_probs)
    else:
        print(f'Gene {Path(f).stem} does not have any high prob. structural homolog')


with open(snakemake.output[0], 'wb') as fh:
    pickle.dump(selected_high, fh)