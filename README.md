# TaxonomicClassification

<!-- markdownlint-disable MD010 -->

This project extends a GraphBin2 metagenomic binning workflow with taxonomic
labels produced by Kraken2. It compares three ways of assigning taxonomy to
reads:

1. Kraken2 with a standard database;
2. Kraken2 with a smaller database;
3. Kraken2 contig labels propagated through GraphBin2 and mapped back to reads.

The workflow is implemented in Snakemake and is designed for CAMI-style
simulated reads, where a gold-standard read-to-taxonomy mapping is available
for evaluation.

## How the workflow works

The default execution follows this data flow:

```text
CAMI reads (.fq/.fq.gz)
    |
    +--> MEGAHIT or MetaSPAdes --> contigs + assembly graph
    |                                  |
    |                                  +--> Kraken2 contig classification
    |                                  |       |
    |                                  |       +--> processor.py --> contig bins
    |                                  |
    |                                  +--> Bowtie2 + Samtools --> abundance + BAM
    |                                                          |
    |                         contig bins + graph + abundance + BAM
    |                                           |
    |                                           +--> GraphBin2
    |                                                   |
    |                                                   +--> reads_binning.py --> read bins
    |
    +--> Kraken2 with standard database --> read classifications
    +--> Kraken2 with mini database     --> read classifications

CAMI truth + all predictions --> generate_binning.py --> AMBER evaluation
```

### Pipeline stages

1. **Assembly**: `megahit` is used by default. Set `assembler: "spades"` in
 `config.yaml` to use MetaSPAdes instead. MEGAHIT output is converted from
 FASTG to GFA before it is passed to GraphBin2. The current rule invokes
 `scripts/fastg2gfa`; the repository also contains `scripts/fastg_to_gfa.py`,
 so this conversion helper should be checked before running a fresh checkout.
2. **Kraken2 classification**: contigs are classified with the mini database;
  reads are classified independently with both the standard and mini
  databases.
3. **Taxonomic bin preparation**: `scripts/processor.py` keeps classified
  records, optionally filters contigs by length and abundance, and walks the
  NCBI taxonomy tree to the configured rank. It writes GraphBin2's contig
  labels plus `bin_taxid_mapping.csv`.
4. **Abundance estimation**: Bowtie2 maps the interleaved reads back to the
  selected contigs. Samtools produces a sorted/indexed BAM and an abundance
  table from `idxstats`.
5. **GraphBin2 propagation**: GraphBin2 uses the assembly graph, contigs,
  abundance table and Kraken-derived contig bins to propagate labels. With
  SPAdes, the `contigs.paths` file is also supplied.
6. **Read assignment**: `scripts/reads_binning.py` converts GraphBin2's
  contig labels to read labels. The workflow currently keeps primary mapped
  reads with MAPQ at least 30.
7. **Evaluation**: `scripts/generate_binning.py` converts the CAMI truth and
  predictions to AMBER `.binning` files. AMBER compares Kraken Standard,
  Kraken Mini and GraphBin results using the configured NCBI taxonomy.

## Repository structure

```text
TaxonomicClassification/
├── Snakefile                 # Snakemake rules and workflow dependencies
├── config.yaml               # Input paths, databases and pipeline settings
├── commands.txt              # Older manual command examples
├── envs/                     # Conda environments used by Snakemake rules
│   ├── assembler.yml         # Assemblers, Kraken2, Bowtie2, Samtools, Python
│   ├── graphbin2.yml         # GraphBin2, Python and pysam
│   ├── amber.yml             # CAMI AMBER
│   ├── opal.yml              # Optional CAMI OPAL environment
│   └── metabat2.yml          # Optional MetaBAT2 environment
├── scripts/
│   ├── processor.py          # Kraken2 output -> GraphBin2 contig bins
│   ├── reads_binning.py      # Contig bins + BAM -> read bins
│   ├── generate_binning.py   # Predictions/truth -> AMBER format
│   ├── evaluate_binning.py   # Direct precision, recall and recovery metrics
│   ├── generate_profile.py   # Profile-generation helper
│   └── fastg_to_gfa.py       # FASTG/GFA conversion helper
└── data/                     # Generated results and local example data
  ├── input_reads/
  ├── megahit_results/
  ├── spades_results/
  ├── kraken_results/
  ├── abundance_results/
  ├── graphBin_results/
  └── benchmark_results/
```

The exact output names are defined by the rules in `Snakefile`. In particular,
the final AMBER directory is `data/benchmark_results/amber_results/` with the
default configuration. The checked-in `data/` directory also contains example
or previously generated results; it is not a replacement for the external
reads and Kraken databases required to run the complete workflow.

## Configuration

Edit `config.yaml` before running the pipeline. The important settings are:

| Key | Meaning | Default |
| --- | --- | --- |
| `assembler` | `megahit` or `spades` | `megahit` |
| `cami_reads` | Interleaved CAMI FASTQ/FASTQ.GZ input | external path |
| `truth_binning` | CAMI gold-standard mapping | external path |
| `kraken_db_mini` | Kraken2 mini database | external path |
| `kraken_db_standard` | Kraken2 standard database | external path |
| `ncbi_taxonomy` | Directory containing `nodes.dmp` | external path |
| `rank` | Taxonomic rank used for binning | `species` |
| `length_threshold` | Minimum contig length | `0` |
| `read_length` | Read length written to AMBER files | `150` |
| `threads` | Maximum threads for multi-threaded rules | `32` |

The paths in the committed configuration point to a sibling `../databases/`
directory. A Kraken2 database directory must contain the database files and
`nodes.dmp`; the taxonomy path used by AMBER must also be compatible with the
database used by the workflow.

## Installation

Install Conda or Micromamba and Snakemake. The remaining rule-specific
dependencies are created automatically from `envs/` when `--use-conda` is
passed.

```bash
conda install -n base -c conda-forge -c bioconda snakemake
```

Check that the external input paths in `config.yaml` exist before starting.
The database checks should be performed according to the files present in your
Kraken2 installation, since database layouts can differ.

## Running the workflow

From the repository root, run the complete workflow with:

```bash
snakemake --use-conda --cores 32
```

On a smaller machine, lower `--cores`; the default configuration asks rules to
use up to 32 threads. To inspect the planned jobs without running them:

```bash
snakemake --use-conda --cores 4 --dry-run
```

To use MetaSPAdes, change the assembler in `config.yaml` and rerun:

```yaml
assembler: "spades"
```

Snakemake will reuse completed outputs when their inputs and parameters have
not changed. Rule-level timing files are written under
`data/benchmark_results/`.

## Important outputs

| Location | Contents |
| --- | --- |
| `data/megahit_results/` or `data/spades_results/` | Contigs and assembly graph |
| `data/kraken_results/` | Kraken2 reports, raw outputs and processed bins |
| `data/abundance_results/` | Bowtie2 index, sorted/indexed BAM and abundance TSV |
| `data/graphBin_results/mini/` | GraphBin2 contig labels and read assignments |
| `data/benchmark_results/binning/` | AMBER input files for truth and predictions |
| `data/benchmark_results/amber_results/` | AMBER comparison reports |

`processor.py` writes a bin-to-taxonomy mapping because GraphBin2 labels are
local bin names, while AMBER needs taxonomy IDs. Do not remove
`bin_taxid_mapping.csv` when inspecting or reusing GraphBin2 read assignments.

## Running helper scripts manually

The scripts can also be used outside Snakemake. For example, to generate a
GraphBin2-compatible contig bin file:

```bash
python3 scripts/processor.py \
 --kraken data/kraken_results/kraken_mini/kraken_output.txt \
 --abundance data/abundance_results/abundance.tsv \
 --contigs data/megahit_results/final_contigs.fa \
 --output data/kraken_results/kraken_mini/output.csv \
 --mapping data/kraken_results/kraken_mini/bin_taxid_mapping.csv \
 --rank species \
 --assembler megahit \
 --db-path /path/to/kraken_db
```

To project those labels onto reads, the BAM must be indexed:

```bash
python3 scripts/reads_binning.py \
 --alignment data/abundance_results/alignment.bam \
 --bin data/kraken_results/kraken_mini/output.csv \
 --o data/graphBin_results/mini/results.csv \
 --minmapq 30
```

For a direct metric report, `evaluate_binning.py` accepts the CAMI truth,
both Kraken read outputs, the GraphBin read output, the bin mapping and
`nodes.dmp`.

## Notes and limitations

- The workflow expects interleaved paired reads through Kraken2, Bowtie2 and
 the assembler.
- `config.yaml` is intentionally environment-specific: paths to reads,
 truth data and databases must be changed on another machine.
- Kraken2 contig classification uses the mini database by default, while
 direct read classification is run against both databases.
- GraphBin2 read assignments are restricted to primary alignments with MAPQ
 30 or higher.
- `envs/opal.yml`, `envs/metabat2.yml` and parts of `commands.txt` document
 optional or earlier experiments; they are not required by the default
 Snakemake target.

## Source

<https://github.com/Pizz0x/TaxonomicClassification>

<!-- markdownlint-enable MD010 -->
