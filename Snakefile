configfile: "config.yaml"

rule all:
    input:
        "data/graphBin_results/graphbin2_output.csv"

rule megahit:
    input:
        r1 = config["reads_1"],
        r2 = config["reads_2"]
    params:
        path = config["megahit_path"]
    output:
        contigs = config["megahit_path"] + "final_contigs.fa",
        input_contigs = config["megahit_path"] + "k141.contigs.fa",
        graph = config["megahit_path"] + "final_graph.gfa"
    shell:
        """
        rm -rf {params.path}tmp

        megahit -1 {input.r1} -2 {input.r2} -o {params.path}tmp

        MAX_K=$(ls {params.path}tmp/intermediate_contigs/k*.final.contigs.fa | sed 's/.*k//;s/.final.contigs.fa//' | sort -n | tail -n 1)

        megahit_toolkit contig2fastg $MAX_K {params.path}tmp/intermediate_contigs/k${{MAX_K}}.contigs.fa > {params.path}tmp/tmp_graph.fastg
        python3 scripts/fastg_to_gfa.py --input {params.path}tmp/tmp_graph.fastg --output {output.graph}

        # megahit_core export-gfa ${{MAX_K_DIR}}/$(basename ${{MAX_K_DIR}}).contigs.fa.gz > {output.graph}

        cp {params.path}tmp/final.contigs.fa {output.contigs}
        cp {params.path}tmp/intermediate_contigs/k${{MAX_K}}.contigs.fa {output.input_contigs}
        """

rule kraken2:
    input:
        contigs = rules.megahit.output.contigs
    output:
        out = f"{config['kraken_path']}/kraken_output.txt",
        report = f"{config['kraken_path']}/kraken_report.txt"
    params:
        db = "~/kraken_db"
    threads: 4
    shell:
        "kraken2 --db {params.db} {input.contigs} > {output.out} --report {output.report} --threads {threads}"

rule abundance:
    input:
        r1 = config["reads_1"],
        r2 = config["reads_2"],
        contigs = rules.megahit.output.contigs
    params:
        index = f"{config['abundance_path']}tmp/contigs_index"
    output:
        alignment = f"{config['abundance_path']}tmp/alignment.bam",
        abundance = f"{config['abundance_path']}abundance.tsv"
    shell:
        """
        bowtie2-build {input.contigs} {params.index}
        bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output.alignment} -
        samtools index {output.alignment}
        samtools idxstats {output.alignment} | grep -v "*" |cut -f 1,3 > {output.abundance}
        """
        


rule preprocessor:
    input:
        kraken = rules.kraken2.output.out,
        abundance = rules.abundance.output.abundance,
        contigs = rules.megahit.output.contigs
    output:
        bins = f"{config['kraken_path']}/output.csv"
    params:
        rank = config["rank"],
        len = config["length_threshold"]
    shell:
        "python3 ./scripts/processor.py --kraken {input.kraken} --abundance {input.abundance} --contigs {input.contigs} --output {output.bins} --rank {params.rank} --len {params.len}"

rule graphbin2:
    input: 
        abundance = rules.abundance.output.abundance,
        contigs = rules.megahit.output.input_contigs,
        graph = rules.megahit.output.graph,
        bins = rules.preprocessor.output.bins
    output:
        final_csv = "data/graphBin_results/graphbin2_output.csv"
    shell:
        "graphbin2 --assembler megahit --graph {input.graph} --contigs {input.contigs} --binned {input.bins} --abundance {input.abundance} --output data/graphBin_results/"