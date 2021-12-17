emapper_pref = config["pref_emapper"]

rule eggnog_mapper:
    input: 
        rules.mmseqs_repr_seqs.output,
    output:
        f"results/eggnog_mapper/{emapper_pref}.emapper.annotations"
    log:
         "logs/emapper.log"
    params:
        db = config["egg_db"],
        prefix = emapper_pref,
        out_dir = "results/eggnog_mapper",
        method = config["method"],
        tax_scope = config["tax_scope"],
        query_cover = config["query_cover"],
        subject_cover = config["subject_cover"],
    threads: config["threads"],
    shell:
        """
        mkdir -p {params.out_dir}
        emapper.py --cpu {threads} -i {input} --output {params.prefix} --output_dir {params.out_dir} --data_dir {params.db} \
        -m {params.method} --sensmode very-sensitive --tax_scope {params.tax_scope} --go_evidence all \
        --target_orthologs all --dbmem --md5 --override --query_cover {params.query_cover} \
        --subject_cover {params.subject_cover} &> {log}
        """

rule select_no_kos:
    input: 
        rules.eggnog_mapper.output
    output:
        f"results/annot_selection/{emapper_pref}_no_kos.list"
    log:
         "logs/select_no_kos.log"
    script:
        "../scripts/egg_select_no_kos.py"
  
rule write_no_ko_seqs:
    input:
        rules.select_no_kos.output,
        rules.write_multi_fasta.output
    output:
        f"results/annot_selection/{emapper_pref}_no_kos.faa"
    log:
         "logs/write_no_ko_seqs.log"
    script:
        "../scripts/write_no_ko_seqs.py"

rule write_subset_seqs:
    input:
        nokos = rules.select_no_kos.output,
        tsv= rules.mmseqs_createtsv.output,
        faa= rules.mmseqs_repr_seqs.output
    output:
        out_dir = directory("results/annot_selection/selected_seqs"),
        seq_list = "results/annot_selection/selected_seqs.list",
        pickle ="results/pickles/rep_memb_dic.pkl"
    log:
         "logs/write_subset_seqs.log"
    params:
        min_n = 1,
        max_n = 100
    threads: 4
    shell:
        """
        mkdir -p {output.out_dir}

        python scripts/write_nokos_seqs_sel.py {input.nokos} {input.tsv} {input.faa} \
        {output.out_dir} {output.seq_list} {threads} {output.pickle} --min_number {params.min_n} \
        --max_number {params.max_n}
        """


checkpoint select_by_size:
    input:
        "results/annot_selection/selected_seqs"
    output:
        directory("results/selected_seqs_by_size")
    params:
        sizes = config["sizes"]
    shell:
        """
        mkdir -p {output[0]} 
        # The second -size is only for testing, remove it for a real run
        find {input[0]} -size {params.sizes} -exec cp {{}} {output[0]} \;
        """