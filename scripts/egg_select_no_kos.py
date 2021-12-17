import sys 

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd


def eggannot_to_df(annot_f):
    '''
    From the eggnog annotation file (`emapper.py`) to a dataframe
    (emapper v2.1.6)
    '''
    import pandas as pd
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

df_nokos = get_nan(snakemake.input[0], 'KEGG_ko')
# This is a list
df_nokos['query'].to_csv(snakemake.output[0], index=False, header=False)
