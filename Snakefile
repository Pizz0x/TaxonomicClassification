configfile: "config.yaml"
#config['kraken_db_path'] = os.path.expanduser(config['kraken_db_path'])
assembler = config['assembler'].lower()

rule all:
    input:
        # config['graphbin_path'] + "graphbin2_output.csv",
        # config['benchmark_path'] + "evaluation_report.txt",
        config['benchmark_path'] + "binning_report.txt"

rule megahit:
    input:
        # r1 = config["reads_1"],
        # r2 = config["reads_2"]
        reads = config['cami_reads']
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
        
        # in case of 2 different reads:  megahit -1 {input.r1} -2 {input.r2} -o {params.path}tmp
        megahit --12 {input.reads} -o {params.path}tmp

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
        # r1 = config['reads_1'],
        # r2 = config['reads_2']
        read = config['cami_reads']
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
        # in case of 2 different reads:  metaspades.py -1 {input.r1} -2 {input.r2} -o {params.path}tmp -t 16
        metaspades.py --12 {input.read} -o {params.path}tmp -t 16

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

rule kraken2_contigs:
    input:
        contigs = lambda wildcards: get_assembler(wildcards)['contigs']
    output:
        out = f"{config['kraken_path']}kraken_mini/kraken_output.txt",
        report = f"{config['kraken_path']}kraken_mini/kraken_report.txt"
    params:
        db = config['kraken_db_mini']
    threads: 4
    conda:
        "envs/assembler.yml"
    shell:
        "kraken2 --db {params.db} {input.contigs} > {output.out} --report {output.report} --threads {threads}"

rule kraken2_reads:
    input:
        # contigs = lambda wildcards: get_assembler(wildcards)['contigs']
        # r1 = config['reads_1'],
        # r2 = config['reads_2']
        reads = config['cami_reads']
    output:
        out = f"{config['kraken_path']}kraken_standard/kraken_reads_output.txt",
        report = f"{config['kraken_path']}kraken_standard/kraken_reads_report.txt"
    params:
        db = config['kraken_db_standard'],
        # Aggiungiamo un check automatico per la compressione
        gzip = "--gzip-compressed" if config['cami_reads'].endswith(".gz") else ""
    threads: 4
    conda:
        "envs/assembler.yml"
    shell:
        """
        kraken2 --db {params.db} \
        {params.gzip} \
        --fastq-input \
        --interleaved \
        --report {output.report} \
        --threads {threads} \
        {input.reads} > {output.out}
        """

rule abundance:
    input:
        # r1 = config["reads_1"],
        # r2 = config["reads_2"],
        reads = config['cami_reads'],
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
        mkdir -p {config[abundance_path]}tmp
        bowtie2-build {input.contigs} {params.index}
        # bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output.alignment} -
        bowtie2 -x {params.index} --interleaved {input.reads} | samtools view -bS - | samtools sort -o {output.alignment} -
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
        kraken = rules.kraken2_contigs.output.out,
        abundance = rules.abundance.output.abundance,
        #db = config['kraken_db_path'],
        contigs = lambda wildcards: get_assembler(wildcards)['contigs']
    output:
        bins = f"{config['kraken_path']}kraken_mini/output.csv",
        mapping = f"{config['kraken_path']}kraken_mini/bin_taxid_mapping.csv"
    params:
        rank = config["rank"],
        len = config["length_threshold"],
        db= config['kraken_db_mini']
    shell:
        "python3 ./scripts/processor.py --kraken {input.kraken} --abundance {input.abundance} --contigs {input.contigs} --output {output.bins} --mapping {output.mapping} --rank {params.rank} --len {params.len} --assembler {assembler} --db-path {params.db}"

# rule preprocessor_reads:
#     input:
#         kraken = rules.kraken2_reads.output.out,
#         db = config['kraken_db_path']
#     output:
#         bins = f"{config['kraken_path']}output_reads.csv"
#     params:
#         rank = config["rank"],
#         len = config["length_threshold"]
#     shell:
#         "python3 ./scripts/processor.py --kraken {input.kraken} --output {output.bins} --rank {params.rank} --len {params.len} --assembler {assembler} --reads --db-path {input.db}"


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
        graphbin_csv = config['graphbin_path'] + "mini/graphbin2_output.csv",
        results = config['graphbin_path'] + "mini/results.csv"
    conda:
        "envs/graphbin2.yml"
    params:
        paths_flag = lambda wildcards, input: f"--paths {input.paths}" if assembler=='spades' else ""
    shell:
        """
        graphbin2 --assembler {assembler} --graph {input.graph} --contigs {input.contigs} --binned {input.bins} --abundance {input.abundance} --output data/graphBin_results/mini/ {params.paths_flag}
        python3 ./scripts/reads_binning.py --alignment {input.alignment} --bin {output.graphbin_csv} --o {output.results} --minmapq 10
        """

# rule benchmark:
#     input:
#         kraken = rules.kraken2_contigs.output.out,
#         graphbin = rules.graphbin2.output.graphbin_csv,
#         abundance = rules.abundance.output.abundance,
#         truth = "./data/mock_truth.tsv",
#         kraken_reads= rules.kraken2_reads.output.out,
#         graphbin_reads=rules.graphbin2.output.results
#     output:
#         report = f"{config['benchmark_path']}evaluation_report_{db_size}.txt"
#     shell:
#         """
#         mkdir -p data/benchmark_results/
#         python3 ./scripts/evaluate_benchmark.py \
#             --kraken {input.kraken} \
#             --graphbin {input.graphbin} \
#             --abundance {input.abundance} \
#             --truth {input.truth} \
#             --kraken_reads {input.kraken_reads} \
#             --graphbin_reads {input.graphbin_reads} > {output.report}
#         """

rule benchmark_binning:
    input:
        truth = config['truth_binning'],
        kraken_standard= rules.kraken2_reads.output.out,
        mapping= rules.preprocessor.output.mapping,
        graphbin_reads= rules.graphbin2.output.results
    output:
        report= config['benchmark_path'] + "binning_report.txt"
    shell:
        """
        python3 scripts/evaluate_binning.py \
            --truth {input.truth} \
            --kraken-reads {input.kraken_standard} \
            --graphbin-reads {input.graphbin_reads} \
            --mapping {input.mapping} > {output.report}
        """