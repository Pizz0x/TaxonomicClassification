# TaxonomicClassification
implementation of an extension for GraphBin2 such that it can accept in input the results obtained from a classifier like Kraken2


Command to execute GraphBin2
```graphbin2 --assembler megahit --graph ./tests/data/ESC_MEGAHIT/final.gfa --contigs ./tests/data/ESC_MEGAHIT/final.contigs.fa --binned ./tests/data/ESC_MEGAHIT/initial_binning_res.csv --abundance ./tests/data/ESC_MEGAHIT/abundance.tsv --output ./tests/data/ESC_MEGAHIT```


Command to execute Kraken2:
```kraken2 --db ~/kraken_db --threads 4 --report ./data/kraken_results/report.txt ./GraphBin2/tests/data/ESC_MEGAHIT/final.contigs.fa > ./data/kraken_results/output.txt```