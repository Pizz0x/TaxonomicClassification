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

        MAX_K=$(ls {params.path}tmp/intermediate_contigs/k*.contigs.fa | sed 's/.*k//;s/.contigs.fa//' | sort -n | tail -n 1)


        megahit_toolkit contig2fastg $MAX_K {params.path}tmp/intermediate_contigs/k${{MAX_K}}.contigs.fa > {output.graph}

        cp {params.path}tmp/final.contigs.fa {output.contigs}
        cp {params.path}tmp/intermediate_contigs/k141.contigs.fa {output.input_contigs}

        rm -rf {params.path}tmp
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

rule preprocessor:
    input:
        kraken = rules.kraken2.output.out,
        abundance = config["abundance_path"],
        contigs = rules.megahit.output.contigs
    output:
        bins = f"{config['kraken_path']}/output.csv"
    params:
        rank = "species",
        len = 600
    shell:
        "python3 ./scripts/processor.py --kraken {input.kraken} --abundance {input.abundance} --contigs {input.contigs} --output {output.bins} --rank {params.rank} --len {params.len}"

rule graphbin2:
    input: 
        abundance = config["abundance_path"],
        contigs = rules.megahit.output.input_contigs,
        graph = rules.megahit.output.graph,
        bins = rules.preprocessor.output.bins
    output:
        final_csv = "data/graphBin_results/graphbin2_output.csv"
    shell:
        "graphbin2 --assembler megahit --graph {input.graph} --contigs {input.contigs} --binned {input.bins} --abundance {input.abundance} --output data/graphBin_results/"