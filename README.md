# Annotating "hypothetical" proteins with the PDB

**See `config/` for configuration information.**

This workflow takes as input a set of protein sequences, clusters them and functionally annotates the clusters' representatives using **Eggnog DB**. Then, it selects those without KO annotations (*"hypothetical proteins"*), aligns them with `hhblits` against **Uniclust30**, and finally, it aligns the resulting MSA against **PDB70**.  

## Testing

Decompress `selected_seqs_by_size.tar.gz` and use that path in the config file (already set).

To see the commands being executed (`-p`) without an actual execution of the workflow, use `-n`. `-r` prints the "reason" for execution of each rule.

```sh
snakemake  --cores 16 -r -p -n
```

`--cores N` specify the max. number of cores used by the whole workflow, so if a rule has set more cores, 
it will use no more than `N`.

**Without the `-n` the workflow will be executed.**

## Results

All the results will be placed inside `/results`. The file **`all_genes_kos.tsv`** presents a list of all the genes which have
one or more KO terms assigned (the rule `propagate_annotations` propagates the annotations from the cluster representatives to 
their members). That file then is used to build a new table, compatible with `ko_mapper.py`, which will produce 3 files:

- `{prefix}_module_completeness.tab`
- `{prefix}_heatmap.pdf`
- `{prefix}_barplot.pdf`


## Others

The `rules.pdf` represents the DAG of this workflow, but it doesn't include the rules related to
`hhblits`, because those rules depend on a `checkpoint` rule. The rules not included in the DAG visualization still are executed, this has to do with the way snakemake works and computes the DAG beforehand, it isn't a bug.

