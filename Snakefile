configfile: "config.yaml"
#config['kraken_db_path'] = os.path.expanduser(config['kraken_db_path'])
assembler = config['assembler'].lower()

rule all:
    input:
        directory(f"{config['benchmark_path']}amber_results")

rule megahit:
    input:
        reads = config['cami_reads']
    params:
        path = config["megahit_path"]
    output:
        contigs = config["megahit_path"] + "final_contigs.fa",
        input_contigs = config["megahit_path"] + "k141.contigs.fa",
        graph_fastg = config["megahit_path"] + "final_graph.fastg"
    benchmark:
        config['benchmark_path'] + "megahit.tsv"
    conda:
        "envs/assembler.yml"
    shell:
        """
        rm -rf {params.path}tmp {output.contigs} {output.input_contigs} {output.graph_fastg}
    
        megahit --12 {input.reads} -o {params.path}tmp

        MAX_K=$(ls {params.path}tmp/intermediate_contigs/k*.final.contigs.fa | sed 's/.*k//;s/.final.contigs.fa//' | sort -n | tail -n 1)

        megahit_toolkit contig2fastg $MAX_K {params.path}tmp/intermediate_contigs/k${{MAX_K}}.contigs.fa > {output.graph_fastg}

        cp {params.path}tmp/final.contigs.fa {output.contigs}
        cp {params.path}tmp/intermediate_contigs/k${{MAX_K}}.contigs.fa {output.input_contigs}

        rm -rf {params.path}tmp
        """

rule preprocess4megahit:
    input:
        graph = rules.megahit.output.graph_fastg
    output:
        out_graph = config["megahit_path"] + "final_graph1.gfa"
    shell:
        """
        scripts/fastg2gfa {input.graph} > {output.out_graph}
        """


rule spades:
    input:
        read = config['cami_reads']
    params:
        path = config["spades_path"]
    output:
        contigs = config['spades_path'] + "contigs.fasta",
        paths = config['spades_path'] + 'contigs.paths',
        graph = config['spades_path'] + 'assembly_graph_with_scaffolds.gfa'
    threads: config['threads']
    benchmark:
        config['benchmark_path'] + "spades.tsv"
    conda:
        "envs/assembler.yml"
    shell:
        """
        rm -rf {params.path}tmp {output.contigs} {output.paths} {output.graph}
        metaspades.py --12 {input.read} -o {params.path}tmp -t {threads} -m 90

        cp {params.path}tmp/contigs.fasta {output.contigs}
        cp {params.path}tmp/contigs.paths {output.paths}
        cp {params.path}tmp/assembly_graph_with_scaffolds.gfa {output.graph}
        rm -rf {params.path}tmp
        """


def get_assembler(wildcards):
    if assembler == 'megahit':
        return {
            "contigs": rules.megahit.output.contigs,
            "graph": rules.preprocess4megahit.output.out_graph,
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
    threads: config['threads']
    benchmark:
        f"{config['benchmark_path']}kraken_mini.tsv"
    conda:
        "envs/assembler.yml"
    shell:
        """
        rm -f {output.out} {output.report}
        kraken2 --db {params.db} {input.contigs} > {output.out} --report {output.report} --threads {threads}
        """

rule kraken2_reads:
    input:
        reads = config['cami_reads']
    output:
        out = f"{config['kraken_path']}kraken_standard/kraken_reads_output.txt",
        report = f"{config['kraken_path']}kraken_standard/kraken_reads_report.txt"
    params:
        db = config['kraken_db_standard'],
        # Aggiungiamo un check automatico per la compressione
        gzip = "--gzip-compressed" if config['cami_reads'].endswith(".gz") else ""
    threads: config['threads']
    benchmark:
        config['benchmark_path'] + "kraken_standard_reads.tsv"
    conda:
        "envs/assembler.yml"
    shell:
        """
        rm -f {output.out} {output.report}
        kraken2 --db {params.db} \
        {params.gzip} \
        --fastq-input \
        --interleaved \
        --report {output.report} \
        --threads {threads} \
        {input.reads} > {output.out}
        """

rule kraken2_mini:
    input:
        reads = config['cami_reads']
    output:
        out = f"{config['kraken_path']}kraken_mini_reads/kraken_reads_output.txt",
        report = f"{config['kraken_path']}kraken_mini_reads/kraken_reads_report.txt"
    params:
        db = config['kraken_db_mini'],
        # Aggiungiamo un check automatico per la compressione
        gzip = "--gzip-compressed" if config['cami_reads'].endswith(".gz") else ""
    threads: config['threads']
    conda:
        "envs/assembler.yml"
    shell:
        """
        rm -f {output.out} {output.report}
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
        reads = config['cami_reads'],
        contigs = lambda wildcards: get_assembler(wildcards)['contigs']
    params:
        index = f"{config['abundance_path']}tmp/contigs_index",
        tmp_dir = f"{config['abundance_path']}tmp"
    output:
        alignment = f"{config['abundance_path']}alignment.bam",
        abundance = f"{config['abundance_path']}abundance.tsv"
    threads: config['threads']
    conda:
        "envs/assembler.yml"
    benchmark:
        f"{config['benchmark_path']}bowtie.tsv"
    shell:
        """
        rm -rf {params.tmp_dir} {output.alignment} {output.alignment}.bai {output.abundance}
        mkdir -p {config[abundance_path]}tmp
        bowtie2-build {input.contigs} {params.index}
        bowtie2 -p {threads} -x {params.index} --interleaved {input.reads} | samtools view -bS - | samtools sort -@ {threads} -o {output.alignment} -
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
        """
        rm -f {output.bins} {output.mapping}
        python3 ./scripts/processor.py --kraken {input.kraken} --abundance {input.abundance} --contigs {input.contigs} --output {output.bins} --mapping {output.mapping} --rank {params.rank} --len {params.len} --assembler {assembler} --db-path {params.db}
        """

rule preprocessor_reads:
    input:
        kraken_reads = rules.kraken2_reads.output.out,
        kraken_mini = rules.kraken2_mini.output.out
    output:
        standard_bins = f"{config['kraken_path']}output_kraken_reads.csv",
        mini_bins = f"{config['kraken_path']}output_kraken_mini.csv"
    params:
        rank = config["rank"],
        len = config["length_threshold"],
        db_standard= config['kraken_db_standard'],
        db_mini=config['kraken_db_mini']
    shell:
        """
        python3 ./scripts/processor.py --kraken {input.kraken_reads} --output {output.standard_bins} --rank {params.rank} --len {params.len} --assembler {assembler} --reads --db-path {params.db_standard}
        python3 ./scripts/processor.py --kraken {input.kraken_mini} --output {output.mini_bins} --rank {params.rank} --len {params.len} --assembler {assembler} --reads --db-path {params.db_mini}
        """


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
    benchmark:
        f"{config['benchmark_path']}graphbin2.tsv"
    params:
        paths_flag = lambda wildcards, input: f"--paths {input.paths}" if assembler=='spades' else "",
        out_dir = config['graphbin_path'] + "mini"
    shell:
        """
        ulimit -n $(ulimit -Hn)
        rm -rf {params.out_dir} {output.results}
        mkdir -p {params.out_dir}
        graphbin2 --assembler {assembler} --graph {input.graph} --contigs {input.contigs} --binned {input.bins} --abundance {input.abundance} --output data/graphBin_results/mini/ {params.paths_flag}
        python3 ./scripts/reads_binning.py --alignment {input.alignment} --bin {output.graphbin_csv} --o {output.results} --minmapq 30
        """

rule amber_evaluation:
    input:
        length = config['read_length'],
        truth = config['truth_binning'],
        kraken_standard = rules.kraken2_reads.output.out,
        kraken_mini = rules.kraken2_mini.output.out,
        mapping = rules.preprocessor.output.mapping,
        graphbin = rules.graphbin2.output.results,
        ncbi = config['ncbi_taxonomy']
    output:
        binning_standard = f"{config['benchmark_path']}binning/kraken_standard.binning",
        binning_mini = f"{config['benchmark_path']}binning/kraken_mini.binning",
        binning_graphbin = f"{config['benchmark_path']}binning/graphbin.binning",
        truth_with_length = f"{config['benchmark_path']}binning/gs_read_mapping_with_length.binning",
        amber_dir = directory(f"{config['benchmark_path']}amber_results")
    conda:
        "envs/amber.yml"
    shell:
        """
        GOLD_SAMPLE=$(zgrep -m 1 '^@SampleID:' {input.truth} | cut -d':' -f2 | tr -d '[:space:]' || echo "Sample_1")
        python3 scripts/generate_binning.py --input {input.truth} --type goldstandard --output {output.truth_with_length} --sample-id "$GOLD_SAMPLE" --length {input.length}

        python3 scripts/generate_binning.py --input {input.kraken_standard} --type kraken --output {output.binning_standard} --sample-id "$GOLD_SAMPLE" --length {input.length}
        python3 scripts/generate_binning.py --input {input.kraken_mini} --type kraken --output {output.binning_mini} --sample-id "$GOLD_SAMPLE" --length {input.length}
        python3 scripts/generate_binning.py --input {input.graphbin} --type graphbin --mapping {input.mapping} --output {output.binning_graphbin} --sample-id "$GOLD_SAMPLE" --length {input.length}

        rm -rf {output.amber_dir}

        amber.py -g {output.truth_with_length} \
             -o {output.amber_dir} \
             -l "Kraken_Standard,Kraken_Mini,GraphBin" \
             --ncbi_dir {input.ncbi} \
             {output.binning_standard} {output.binning_mini} {output.binning_graphbin}
        """