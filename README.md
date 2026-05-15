# TaxonomicClassification

Implementation of an extension for GraphBin2 such that it can accept in input the results obtained from a classifier like Kraken2

### Setting Up the Repository

To obtain the envirmonment required to run the code you need to execute: \
```conda env create -f environment.yml```\
Which will create a new conda environment called bio-dev that has all the necessary dependencies.

Then you can run the pipeline using the snakefile, using the command:\
```snakemake --cores 4```

Possible Next Steps:

- compare kraken with other technologies like CLARK
- dobbiamo trovare il modo per mappare le reads ai contigs dall'assemblatore (forse bowtie dato che usa reads e contigs)
