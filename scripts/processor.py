<<<<<<< HEAD
import csv
import os
import sys
import argparse

"""
This script is used to process Kraken2 output and obtain GraphBin2 binning input.

Required input formats:
    - Kraken2 output: tab-delimited lines in standard Kraken2 format
      the columns we are interested in are the first three: status, sequence ID and taxonomy ID
      only the instances with status = 'C' are keeped
    - FASTA Contigs: contig ID and sequence length are taken from the header line (after '>')
    - Abundance TSV: file with contig ID in the first column.

The script validates required files before processing and rejects invalid rank values early.
"""



rank_value = {
    'superkingdom': 10,
    'kingdom': 9,
    'phylum': 8,
    'class': 7,
    'order': 6,
    'family': 5,
    'genus': 4,
    'species': 3,
    'no rank': 0
}
VALID_RANKS = set(rank_value.keys())


def file_exists(path, desc):
    if not path:
        raise ValueError(f"{desc} path is required")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{desc} not found: {path}")
    return path


def validate_args(args, NCBI, parser):
    try:
        file_exists(args.kraken, 'Kraken output')
        file_exists(NCBI, 'NCBI taxonomy nodes.dmp')
        if not args.reads:
            file_exists(args.abundance, 'Abundance file')
            file_exists(args.contigs, 'Contigs FASTA file')
    except (ValueError, FileNotFoundError) as exc:
        parser.error(str(exc))

    if args.rank not in VALID_RANKS:
        parser.error(f"Invalid rank '{args.rank}'. Valid ranks are: {', '.join(VALID_RANKS)}")
=======
import csv, os, argparse

tax_set = set()
NCBI_taxonomy = os.path.expanduser('~/kraken_db/nodes.dmp')
abundance_path = './GraphBin2/tests/data/ESC_MEGAHIT/abundance.tsv'
contigs_path = './GraphBin2/tests/data/ESC_MEGAHIT/final.contigs.fa'
rank_value = {
    'superkingdom' : 10,
    'kingdom' : 9,
    'phylum' : 8,
    'class' : 7,
    'order' : 6,
    'family' : 5,
    'genus' : 4,
    'species' : 3,
    'no rank' : 0
}

# get kraken output and make it a pair of contig id and taxonomy id (only when the instance is correctly classified)
def process_kraken(kraken_output, abundance_path, contigs_path, min_length):
    raw_dict = {}
    abundance_set = load_abundance(abundance_path)
    contigs_set = load_contigs(contigs_path, min_length)
    with open(kraken_output, 'r') as i:
        for line in i:
            features = line.strip().split('\t')
            status = features[0] # classified or not
            seqid = features[1] # identifier of the contig
            taxid = int(features[2]) # identifier of the taxonomy classification (it can be species, genus, family, order, class, phylum, kingdom or domain), we need to see if selecting a certain generalization improve results
            if status == 'C' and seqid in abundance_set and seqid in contigs_set: # only consider classified contigs
                if taxid not in tax_set:
                    tax_set.add(taxid) # add new taxonomy id to the set
                raw_dict[seqid] = taxid # pairs: contig identifier, taxonomy identifier (raw since we still have to generalize)
    return raw_dict
>>>>>>> df16b5496842f1b4ed9efec35bddacaf51103357

# we use the NCBI taxonomy database to get the generalization of the taxonomy classification, which is stored in the nodes.dmp file inside of the miniDB of Kraken
def load_NCBI(NCBI_taxonomy):
    taxonomy_dict = {}
    with open(NCBI_taxonomy, 'r') as i:
        for line in i:
            features = line.strip().split('\t|\t')
            current = int(features[0])
            parent = int(features[1])
            rank = features[2].lower().strip()
            taxonomy_dict[current] = (parent, rank)
    return taxonomy_dict

def load_abundance(abundance_path):
    abundance_set = set()
    with open(abundance_path, 'r') as i:
        for instance in i:
            features = instance.split('\t')
            abundance_set.add(features[0])
    return abundance_set

# we take all the contigs over a certain length
<<<<<<< HEAD
def load_contigs(contigs_path, min_length=0):
    contigs_set = set()

    with open(contigs_path, 'r') as i:
        for line in i:
            line = line.rstrip('\n')
            if line.startswith('>'):
                node = line[1:].split()[0] # get the identifier of the contig (without > and without the rest of the description)
                try:
                    len = int(node.split('_')[3])
                except IndexError:
                    len = 0
=======
def load_contigs(contigs_path, min_length = 0):
    contigs_set = set()
    with open(contigs_path, 'r') as i:
        for line in i:
            if line.startswith('>'):
                node = line[1:].split()[0] # get the identifier of the contig (without > and without the rest of the description)
                len = int(node.split('_')[3])
>>>>>>> df16b5496842f1b4ed9efec35bddacaf51103357
                if len >= min_length:
                    contigs_set.add(node) # no character >
    return contigs_set

# given a taxonomy identifier, climb the tree until the wanted generalization is obtained
def get_parent(taxid, taxonomy_dict, rank):
<<<<<<< HEAD
    if taxid not in taxonomy_dict:
        return taxid
    current_id = taxid
    
=======
    current_id = taxid
>>>>>>> df16b5496842f1b4ed9efec35bddacaf51103357
    parent_id = taxonomy_dict[current_id][0]
    current_rank = taxonomy_dict[current_id][1]
    current_rank = 'no rank' if current_rank not in rank_value.keys() else current_rank
    while rank_value[current_rank] < rank_value[rank] and current_id != parent_id:
        current_id = parent_id
<<<<<<< HEAD
        if current_id not in taxonomy_dict:
            break
=======
>>>>>>> df16b5496842f1b4ed9efec35bddacaf51103357
        parent_id = taxonomy_dict[current_id][0]
        current_rank = taxonomy_dict[current_id][1]
        current_rank = 'no rank' if current_rank not in rank_value.keys() else current_rank
    return current_id
    

<<<<<<< HEAD
def bin_format(bin_id, assembler):
    if assembler == "megahit":
        label = f"res.00{bin_id}]" if bin_id < 10 else f"res.0{bin_id}" if bin_id < 100 else f"res.{bin_id}"
    elif assembler == "spades":
        label = f"bin_{bin_id}"
    else:
        label = f"bin.{bin_id}"
    return label



# get kraken output and make it a pair of contig id and taxonomy id (only when the instance is correctly classified)
def process_kraken(kraken_bins, abundance, contigs, taxonomy_dict, rank, output_file, assembler, min_length, is_reads=False):

    abundance_set = load_abundance(abundance) if not is_reads else None
    contigs_set = load_contigs(contigs, min_length) if not is_reads else None

    bins = {}
    bin_counter = 1

    parent_lookup = {}

    with open(kraken_bins, 'r') as i, open(output_file, 'w') as o:
        writer = csv.writer(o)

        for line_num, line in enumerate(i, 1):

            features = line.strip().split('\t')
            if len(features) < 3:
                continue
            status = features[0]
            seqid = features[1] # identifier of the contig
            try:
                taxid = int(features[2])  # identifier of the taxonomy classification
            except ValueError:
                print(f"Warning: invalid taxonomy ID on line {line_num}: {features[2]}", file=sys.stderr)
                continue

            if status == 'C' and (is_reads or (seqid in abundance_set and seqid in contigs_set)): # only consider correctly classified instances and filter by abundance and contig length if not using read-level output
                if taxid in parent_lookup:
                    gen_taxid = parent_lookup[taxid]
                else:
                    gen_taxid = get_parent(taxid, taxonomy_dict, rank)
                    parent_lookup[taxid] = gen_taxid
                
                if gen_taxid not in bins:
                    bins[gen_taxid] = bin_counter
                    bin_counter += 1
                
                bin_number = bins[gen_taxid]
                label = bin_format(bin_number, assembler)

                writer.writerow([seqid, label])



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Kraken2 output and create GraphBin2 binned input.")
    parser.add_argument("--rank", default="genus", choices=VALID_RANKS,
                        help="Taxonomic rank to propagate (one of: %(choices)s)")
    parser.add_argument("--len", type=int, default=0,
                        help="Minimum contig length to include when using contig-based Kraken output.")
    parser.add_argument("--kraken", default='./data/kraken_results/output.txt',
                        help="Kraken2 output file with status, sequence ID, taxonomy ID, ...")
    parser.add_argument("--output", default='./data/kraken_results/output.csv',
                        help="CSV file to write contig-to-bin labels.")
    parser.add_argument("--abundance", help="Abundance TSV file with contig IDs in the first column.")
    parser.add_argument("--contigs", help="Contigs FASTA file used to filter contigs by length.")
    parser.add_argument("--assembler", default="megahit",
                        help="Assembler type used for label formatting (megahit or spades).")
    parser.add_argument("--reads", action='store_true',
                        help="Use Kraken2 read-level output instead of contig-level output.")
    parser.add_argument("--db-path", default="~/kraken_db", help="Path to your Kraken database node files.")
    args = parser.parse_args()

    NCBI_taxonomy = os.path.expanduser(f"{args.db_path}/nodes.dmp")
    validate_args(args, NCBI_taxonomy, parser)
    taxonomy_tree = load_NCBI(NCBI_taxonomy)
    process_kraken(args.kraken, args.abundance, args.contigs, taxonomy_tree, args.rank, args.output, args.assembler, args.len, args.reads)
=======
# given the pair of contig_id and taxonomy_id , the dictionary that map each required tax_id to his parent and write each pair on the csv file
# this is the funciton that creates the input for GraphBin2
def process_raw(raw_dict, taxonomy_dict, rank, output_file, assembler):
    bins = {}
    count = 1
    for taxid in raw_dict.values():
        gen_taxid = get_parent(taxid, taxonomy_dict, rank)
        if gen_taxid not in bins:
            bins[gen_taxid] = count
            count += 1

    with open(output_file, 'w') as o:
        writer = csv.writer(o)
        count = 1
        for seqid, taxid in raw_dict.items():
            gen_taxid = get_parent(taxid, taxonomy_dict, rank)
            if assembler == "megahit":
                label = f"res.00{bins[gen_taxid]}" if bins[gen_taxid] < 10 else f"res.0{bins[gen_taxid]}" if bins[gen_taxid] < 100 else f"res.{bins[gen_taxid]}"
            elif assembler == "spades":
                label = f"bin_{bins[gen_taxid]}"
            #label = f"bin.{bins[gen_taxid]}.fa"
            #contig = f"contig_{count}"
            #count +=1
            writer.writerow([seqid, label])
        

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--rank", default="genus")
    parser.add_argument("--len", default=0)
    parser.add_argument("--kraken", default='./data/kraken_results/output.txt')
    parser.add_argument("--output", default='./data/kraken_results/output.csv')
    parser.add_argument("--abundance", default=abundance_path)
    parser.add_argument("--contigs", default=contigs_path)
    parser.add_argument("--assembler", default="megahit")
    args = parser.parse_args()

    raw_dict = process_kraken(args.kraken, args.abundance, args.contigs, int(args.len))
    taxonomy_dict = load_NCBI(NCBI_taxonomy)
    process_raw(raw_dict, taxonomy_dict, args.rank, args.output, args.assembler)
>>>>>>> df16b5496842f1b4ed9efec35bddacaf51103357


#graphbin2 --assembler megahit --graph ./data/megahit_results/final_graph.gfa --contigs ./data/megahit_results/final_contigs.fa --binned ./data/kraken_results/output.csv --abundance ./data/abundance_results/abundance.tsv --output data/graphBin_results/