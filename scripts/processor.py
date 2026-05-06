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
def load_contigs(contigs_path, min_length = 0):
    contigs_set = set()
    with open(contigs_path, 'r') as i:
        for line in i:
            if line.startswith('>'):
                features = line.split()
                len = int(features[3].split('=')[1])
                if len >= min_length:
                    contigs_set.add(features[0][1:]) # no character >
    return contigs_set

# given a taxonomy identifier, climb the tree until the wanted generalization is obtained
def get_parent(taxid, taxonomy_dict, rank):
    current_id = taxid
    parent_id = taxonomy_dict[current_id][0]
    current_rank = taxonomy_dict[current_id][1]
    current_rank = 'no rank' if current_rank not in rank_value.keys() else current_rank
    while rank_value[current_rank] < rank_value[rank] and current_id != parent_id:
        current_id = parent_id
        parent_id = taxonomy_dict[current_id][0]
        current_rank = taxonomy_dict[current_id][1]
        current_rank = 'no rank' if current_rank not in rank_value.keys() else current_rank
    return current_id
    

# given the pair of contig_id and taxonomy_id , the dictionary that map each required tax_id to his parent and write each pair on the csv file
# this is the funciton that creates the input for GraphBin2
def process_raw(raw_dict, taxonomy_dict, rank, output_file):
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
            label = f"res.00{bins[gen_taxid]}" if bins[gen_taxid] < 10 else f"res.0{bins[gen_taxid]}" if bins[gen_taxid] < 100 else f"res.{bins[gen_taxid]}"
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
    args = parser.parse_args()

    raw_dict = process_kraken(args.kraken, args.abundance, args.contigs, int(args.len))
    taxonomy_dict = load_NCBI(NCBI_taxonomy)
    process_raw(raw_dict, taxonomy_dict, args.rank, args.output)


#graphbin2 --assembler megahit --graph ./data/megahit_results/final_graph.gfa --contigs ./data/megahit_results/final_contigs.fa --binned ./data/kraken_results/output.csv --abundance ./data/abundance_results/abundance.tsv --output data/graphBin_results/