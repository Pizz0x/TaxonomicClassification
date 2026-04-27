# TaxonomicClassification
implementation of an extension for GraphBin2 such that it can accept in input the results obtained from a classifier like Kraken2


Command to execute GraphBin2
```graphbin2 --assembler megahit --graph ./GraphBin2/tests/data/ESC_MEGAHIT/final.gfa --contigs ./GraphBin2/tests/data/ESC_MEGAHIT/final.contigs.fa --binned ./data/kraken_results/output.csv --abundance ./GraphBin2/tests/data/ESC_MEGAHIT/abundance.tsv --output ./GraphBin2/tests/data/ESC_MEGAHIT```

```graphbin2 --assembler flye --graph ./GraphBin2/tests/data/1Y3B_Flye/assembly_graph.gfa --contigs ./GraphBin2/tests/data/1Y3B_Flye/assembly.fasta --paths ./GraphBin2/tests/data/1Y3B_Flye/assembly_info.txt --binned ./data/kraken_results/output.csv --abundance ./GraphBin2/tests/data/1Y3B_Flye/abundance.tsv --output ./data```

Command to execute Kraken2:
```kraken2 --db ~/kraken_db --threads 4 --report ./data/kraken_results/report.txt ./GraphBin2/tests/data/ESC_MEGAHIT/final.contigs.fa > ./data/kraken_results/output.txt```