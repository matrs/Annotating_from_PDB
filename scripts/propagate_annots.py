import pandas as pd
from pathlib import Path
import pickle
from collections import defaultdict

def eggannot_to_df(annot_f):
    '''
    From the eggnog annotation file (`emapper.py`) to a dataframe
    (emapper v2.1.6)
    '''
# ParserWarning: Falling back to the 'python' engine because the 'c' engine does 
# not support skipfooter; you can avoid this warning by specifying engine='python'.
    annot_df =  pd.read_csv(annot_f, sep='\t', skiprows=4, skipfooter=3,  na_values='-', engine='python')
    annot_df.rename({'#query':'query'}, axis=1, inplace=True)
    return annot_df

def get_nan(annot_file, colname):
    '''
    Get the NAs values in the column colname
    annot_file: egnogg anotation file
    colname: emapper annotation column name, e.g. 'COG_category', 'KEGG_ko', etc
    '''
    
    df = eggannot_to_df(annot_file)
    mask = df[colname].isna()
    new_df = df[mask]
    return new_df

def get_not_nan(annot_file, colname):
    '''
    Get the "not NAs" values in the column colname
    annot_file: egnogg anotation file
    colname: emapper annotation column name, e.g. 'COG_category', 'KEGG_ko', etc
    '''
    
    df = eggannot_to_df(annot_file)
    mask = df[colname].notna()
    new_df = df[mask]
    return new_df

with open(snakemake.input.pickle, 'rb') as fh:
    map_dic_hyp = pickle.load(fh)

print(repr(snakemake.input), snakemake.input.annot[0], type(snakemake.input.pickle))
ko_annot_df = get_not_nan(snakemake.input[0], 'KEGG_ko')
print("Dimension of annotated KOs DF:", ko_annot_df.shape)

all_genes_annot = defaultdict()
for i, row in ko_annot_df.iterrows():
    for gene in map_dic_hyp[row['query']]:
        all_genes_annot[gene] = row['KEGG_ko']
print(len(all_genes_annot))


print(snakemake.input.pdb_annot)
pdb_annots =  pd.read_csv(snakemake.input.pdb_annot, sep='\t')
mask = pdb_annots['KO'].notna()
pdb_annots = pdb_annots[mask]

for i, row in pdb_annots.iterrows():
    for gene in map_dic_hyp[row['gene']]:
        all_genes_annot[gene] = row['KO']
len(all_genes_annot)

with open(snakemake.output[0], 'a') as fh:
    for gene, kos in all_genes_annot.items():
        print(gene, kos, file=fh, sep='\t')
