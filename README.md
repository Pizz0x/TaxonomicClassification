# TaxonomicClassification

implementation of an extension for GraphBin2 such that it can accept in input the results obtained from a classifier like Kraken2

Command to execute GraphBin2
```graphbin2 --assembler megahit --graph ./data/megahit_results/final_graph.gfa --contigs ./data/megahit_results/k141.contigs.fa --binned ./data/kraken_results/output.csv --abundance ./GraphBin2/tests/data/ESC_MEGAHIT/abundance.tsv --output ./data/graphBin_results```

```graphbin2 --assembler flye --graph ./GraphBin2/tests/data/1Y3B_Flye/assembly_graph.gfa --contigs ./GraphBin2/tests/data/1Y3B_Flye/assembly.fasta --paths ./GraphBin2/tests/data/1Y3B_Flye/assembly_info.txt --binned ./data/kraken_results/output.csv --abundance ./GraphBin2/tests/data/1Y3B_Flye/abundance.tsv --output ./data```

Command to execute Kraken2:
```kraken2 --db ~/kraken_db --threads 4 --report ./data/kraken_results/report.txt ./GraphBin2/tests/data/ESC_MEGAHIT/final.contigs.fa > ./data/kraken_results/output.txt```

another thing that we should try is to compare kraken with other technologies like CLARK

rm -rf ./data/megahit_results_tmp

megahit -1 ./data/input_reads/SRR5058924_1.fastq.gz -2 ./data/input_reads/SRR5058924_2.fastq.gz -o ./data/megahit_results_tmp

MAX_K=$(ls ./data/megahit_results_tmp/intermediate_contigs/k*.contigs.fa | sed 's/.*k//;s/.contigs.fa//' | sort -n | tail -n 1)

megahit_toolkit contig2fastg 141 ./data/megahit_results/final_contigs.fa > ./data/megahit_results/final_graph.gfa

cp ./data/megahit_results_tmp/final.contigs.fa ./data/megahit_results/final.contigs.fa
cp ./data/megahit_results_tmp/final_graph.gfa ./data/megahit_results/final_graph.gfa

rm -rf ./data/megahit_results_tmp
