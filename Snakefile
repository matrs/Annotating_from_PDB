from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("6.12.3")


##### setup report #####
configfile: "config/config.yaml"

include: "rules/get_prot_seqs.smk"
include: "rules/clustering.smk"
include: "rules/eggmapper.smk"
include: "rules/hhblits_to_ko.smk"


prefix = config['pref_mmseqs']
emapper_pref = config["pref_emapper"]
ko_pref = config["ko_prefix"]

rule all:
    input:
         "results/clustering/seqs_counts.txt",
         f"results/clustering/{prefix}_repr.tsv",
         f"results/annot_selection/{emapper_pref}_no_kos.faa",
         "results/annot_selection/selected_seqs.list",
         "results/akk_not_anot_pdbs_kos.tsv",
         "results/all_genes_kos.tsv",
         f"results/{ko_pref}_module_completeness.tab",


