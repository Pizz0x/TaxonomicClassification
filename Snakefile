configfile: "config.yaml"
assembler = config['assembler'].lower()

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
    conda:
        "envs/assembler.yml"
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

        rm -rf {params.path}tmp
        """


rule spades:
    input:
        r1 = config['reads_1'],
        r2 = config['reads_2']
    params:
        path = config["spades_path"]
    output:
        contigs = config['spades_path'] + "contigs.fasta",
        paths = config['spades_path'] + 'contigs.paths',
        graph = config['spades_path'] + 'assembly_graph_with_scaffolds.gfa'
    conda:
        "envs/assembler.yml"
    shell:
        """
        rm -rf {params.path}tmp
        metaspades.py -k 21,33,55,77,99,127 -1 {input.r1} -2 {input.r2} -o {params.path}tmp -t 16

        cp {params.path}tmp/contigs.fasta {output.contigs}
        cp {params.path}tmp/contigs.paths {output.paths}
        cp {params.path}tmp/assembly_graph_with_scaffolds.gfa {output.graph}
        rm -rf {params.path}tmp
        """


def get_assembler(wildcards):
    if assembler == 'megahit':
        return {
            "contigs": rules.megahit.output.contigs,
            "graph": rules.megahit.output.graph,
            "input_contigs": rules.megahit.output.input_contigs,
            "paths": []
        }
    elif assembler == 'spades':
        return {
            "contigs": rules.spades.output.contigs,
            "graph": rules.spades.output.graph,
            "input_contigs": rules.spades.output.contigs,
            "paths": rules.spades.output.paths
        }
    else:
        raise ValueError(f"Unsupported assembler: {assembler}")


rule kraken2:
    input:
        contigs = lambda wildcards: get_assembler(wildcards)['contigs']
    output:
        out = f"{config['kraken_path']}kraken_output.txt",
        report = f"{config['kraken_path']}kraken_report.txt"
    params:
        db = "~/kraken_db"
    threads: 4
    conda:
        "envs/assembler.yml"
    shell:
        "kraken2 --db {params.db} {input.contigs} > {output.out} --report {output.report} --threads {threads}"


rule abundance:
    input:
        r1 = config["reads_1"],
        r2 = config["reads_2"],
        contigs = lambda wildcards: get_assembler(wildcards)['contigs']
    params:
        index = f"{config['abundance_path']}tmp/contigs_index"
    output:
        alignment = f"{config['abundance_path']}alignment.bam",
        abundance = f"{config['abundance_path']}abundance.tsv"
    conda:
        "envs/assembler.yml"
    shell:
        """
        bowtie2-build {input.contigs} {params.index}
        bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output.alignment} -
        samtools index {output.alignment}
        samtools idxstats {output.alignment} | grep -v "*" |cut -f 1,3 > {output.abundance}
        """


# rule metabat2:
#     input:
#         contigs = lambda wildcards: get_assembler(wildcards)['contigs'],
#         alignment = rules.abundance.output.alignment
#     output:
#         bins = f"{config['metabat_path']}bins/",
#         depth = f"{config['metabat_path']}depth.txt"
#     conda:
#         "envs/metabat2.yml"
#     shell:
#         """
#         jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.alignment} 
#         metabat2 -i {input.contigs} -a {output.depth} -o {output.bins}bin -m 1500
#         """

rule preprocessor:
    input:
        kraken = rules.kraken2.output.out,
        abundance = rules.abundance.output.abundance,
        contigs = lambda wildcards: get_assembler(wildcards)['contigs']
    output:
        bins = f"{config['kraken_path']}output.csv"
    params:
        rank = config["rank"],
        len = config["length_threshold"]
    shell:
        "python3 ./scripts/processor.py --kraken {input.kraken} --abundance {input.abundance} --contigs {input.contigs} --output {output.bins} --rank {params.rank} --len {params.len} --assembler {assembler}"


rule graphbin2:
    input: 
        abundance = rules.abundance.output.abundance,
        alignment = rules.abundance.output.alignment,
        contigs = lambda wildcards: get_assembler(wildcards)["input_contigs"],
        graph = lambda wildcards: get_assembler(wildcards)["graph"],
        paths = lambda wildcards: get_assembler(wildcards)["paths"],
        # an alternative is to use metabat2 bins as input
        bins = rules.preprocessor.output.bins
    output:
        graphbin_csv = config['graphbin_path'] + "graphbin2_output.csv",
        results = config['graphbin_path'] + "results.csv"
    conda:
        "envs/graphbin2.yml"
    params:
        paths_flag = lambda wildcards, input: f"--paths {input.paths}" if assembler=='spades' else ""
    shell:
        """
        graphbin2 --assembler {assembler} --graph {input.graph} --contigs {input.contigs} --binned {input.bins} --abundance {input.abundance} --output data/graphBin_results/ {params.paths_flag}
        python3 ./scripts/reads_binning.py --alignment {input.alignment} --bin {output.graphbin_csv} --o {output.results}
        """