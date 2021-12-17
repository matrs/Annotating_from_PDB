from os.path import join

def get_selected_files(wildcards):
    ck_output = checkpoints.select_by_size.get(**wildcards).output[0]
    GENES,  = glob_wildcards(join(ck_output, "{gene}.fasta"))
    return expand("results/hhblits_pdb/{GENE}.hhr", GENE=GENES)


rule hhblits_msa:
    input:
        "results/selected_seqs_by_size/{gene}.fasta"
    output:
        "results/hhblits_msa/{gene}.a3m"
    log:
        "logs/hhblits/{gene}_msa.log" 
    params:
        db=config["msa_db"],
        n_iter= 1
    threads: 4
    shell:
        """
        mkdir -p results/hhblits_msa
        hhblits -i {input} -d {params.db} -oa3m {output} -cpu {threads}  -n {params.n_iter} \
        &> {log}
        """

rule hhblits_pdb:
    input:
        "results/hhblits_msa/{gene}.a3m"
    output:
        "results/hhblits_pdb/{gene}.hhr",
    log:
        "logs/hhblits/{gene}_pdb.log" 
    params:
        db=config["pdb_db"],
        n_iter= 1
    threads: 4
    shell:
        """
        mkdir -p results/hhblits_pdb
        hhblits -i {input[0]} -d {params.db} -o {output[0]} -cpu {threads}  -n {params.n_iter} \
        &> {log}
        """


rule select_pdbs:
    input:
        get_selected_files
    output:
        "results/pickles/high_prob.pkl"
    log:
        "logs/select_pdbs.log" 
    params:
        # out_dir= "results/",
        prob= 70
    script:
        "../scripts/hhblits_select_pdbs_sm.py"


rule pdb_to_ko:
    input:
        "results/pickles/high_prob.pkl"
    output:
        "results/akk_not_anot_pdbs_kos.tsv"
    log:
        "logs/pdb_to_ko.log"
    params:
        out_dir="results/"
    shell:
        """
        python  scripts/pdb_to_kegg.py {input} {output} &> {log}
        """

rule propagate_annotations:
    input:
        annot = rules.eggnog_mapper.output,
        pickle = rules.write_subset_seqs.output.pickle,
        pdb_annot = "results/akk_not_anot_pdbs_kos.tsv"
    output:
        "results/all_genes_kos.tsv",
    log:
        "logs/propagate_annotations.log"
    params:
        out_dir="results/"
    script:
        "../scripts/propagate_annots.py"

rule make_ko_list:
    input:
        rules.propagate_annotations.output,
    output:
        "results/kos.list",
    log:
        "logs/propagate_annotations.log"
    params:
        out_dir="results/"
    script:
        "../scripts/make_ko_list.py"



py_komapper= config["ko_mapper_path"]
ko_pref = config["ko_prefix"]

rule ko_completness:
    input:
        rules.make_ko_list.output,
    output:
        f"results/{ko_pref}_module_completeness.tab",
        f"results/{ko_pref}_heatmap.pdf",
        f"results/{ko_pref}_barplot.pdf",
    log:
        "logs/propagate_annotations.log"
    params:
        out_dir="results/", 
        prefix= ko_pref,
        ko_mapper_path = config["ko_mapper_path"]
    shell:
        "python {py_komapper} -i {input} -p results/{params.prefix} --cluster rows"

