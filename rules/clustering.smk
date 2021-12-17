rule mmseqs_db:
    input: 
        rules.write_multi_fasta.output
    output: 
        seq_db=directory("results/clustering/mmseq_seqdb/"),
    log: "logs/mmseqs_db.log"
    params: 
        prefix=config["pref_mmseqs"],
    shell:
        """
        mkdir -p {output.seq_db}
        mmseqs createdb {input[0]} {output.seq_db}/{params.prefix} &> {log}
        """

rule mmseqs_clustering:
    input: 
        rules.mmseqs_db.output.seq_db,
    output: 
        cluster_db=directory("results/clustering/mmseq_clusterdb/")
    log:
        "logs/mmseqs_cluster.log"
    params: 
        prefix=config["pref_mmseqs"],
        cov=config["coverage"],
        identity=config["identity"]
    threads: 4
    shell:
        """
        mkdir -p {output.cluster_db}
        mmseqs cluster {input[0]}/{params.prefix} {output.cluster_db}/{params.prefix} \
        tmp  --cov-mode 0 -c {params.cov} --min-seq-id {params.identity} \
        --threads {threads} &>> {log}

        # remove the tmp directory
        rm -rf tmp 
        """

rule mmseqs_createsubdb:
    input: 
        rules.mmseqs_clustering.output,
        rules.mmseqs_db.output
    output:
        subdb=directory("results/clustering/mmseq_subdb/"),
    log:
         "logs/mmseqs_subdb.log"
    params: 
        prefix=config["pref_mmseqs"],
    shell:
        """
        mkdir -p {output.subdb}
        mmseqs createsubdb {input[0]}/{params.prefix} {input[1]}/{params.prefix} \
        {output.subdb}/{params.prefix} &> {log}
        """

rule mmseqs_repr_seqs:
    input: 
        rules.mmseqs_createsubdb.output.subdb,
    output:
        f"results/clustering/{prefix}_repr.faa"
    log:
         "logs/mmseqs_repr_seq.log"
    params: 
        prefix=config["pref_mmseqs"],
    shell:
        """
        mmseqs convert2fasta {input[0]}/{params.prefix} {output[0]}   &> {log}
        """


rule mmseqs_createtsv:
    input:
        rules.mmseqs_db.output,
        rules.mmseqs_clustering.output
    output:
       f"results/clustering/{prefix}_repr.tsv"
    log:
         "logs/mmseqs_repr_seq.log"
    params: 
        prefix=config["pref_mmseqs"],
    shell:
        """
        mmseqs createtsv {input[0]}/{params.prefix} {input[0]}/{params.prefix} \
        {input[1]}/{params.prefix} {output} &> {log}
        """


rule count_seqs:
    input: 
        rules.write_multi_fasta.output,
        rules.mmseqs_repr_seqs.output,
    output:
        "results/clustering/seqs_counts.txt"
    shell:
        """
        echo "Number of sequences in the original fasta and the clustered fasta:" >  {output}
        egrep -c '^>' {input[0]} >> {output}
        egrep -c '^>' {input[1]} >> {output}
        """