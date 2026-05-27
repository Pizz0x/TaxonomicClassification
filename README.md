# TaxonomicClassification

Implementation of an extension for GraphBin2 such that it can accept in input the results obtained from a classifier like Kraken2

### Setting Up the Repository

To obtain the envirmonment required to run the code you need to execute: \
```conda env create -f environment.yml```\
Which will create a new conda environment called bio-dev that has all the necessary dependencies.

Then you can run the pipeline using the snakefile, using the command:\
```snakemake --cores 4 --use-conda```

# Next Steps

- implement a processor for the reads classified by kraken2
- check if the fact that we use contigs to obtain the bins and then reconvert it into reads thanks to the use of bowtie is better than directly use the reads for kraken
- check if this works also for kraken with a larger database
- compare kraken with other technologies like CLARK
