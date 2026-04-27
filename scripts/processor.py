import csv, os, argparse

tax_set = set()
NCBI_taxonomy = os.path.expanduser('~/kraken_db/nodes.dmp')
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

def process_kraken(kraken_output):
    raw_dict = {}
    with open(kraken_output, 'r') as i:
        for line in i:
            features = line.strip().split('\t')
            status = features[0] # classified or not
            seqid = features[1] # identifier of the contig
            taxid = int(features[2]) # identifier of the taxonomy classification (it can be species, genus, family, order, class, phylum, kingdom or domain), we need to see if selecting a certain generalization improve results

            if status == 'C': # only consider classified contigs
                if taxid not in tax_set:
                    tax_set.add(taxid) # add new taxonomy id to the set
                raw_dict[seqid] = taxid # pairs: contig identifier, taxonomy identifier (raw since we still have to generalize)
    return raw_dict

def load_NCBI(NCBI_taxonomy):
    # we use the NCBI taxonomy database to get the generalization of the taxonomy classification, which is stored in the nodes.dmp file inside of the miniDB of Kraken
    taxonomy_dict = {}
    with open(NCBI_taxonomy, 'r') as i:
        for line in i:
            features = line.strip().split('\t|\t')
            current = int(features[0])
            parent = int(features[1])
            rank = features[2].lower().strip()
            taxonomy_dict[current] = (parent, rank)
    return taxonomy_dict

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
    



def process_raw(raw_dict, taxonomy_dict, rank, output_file):
    with open(output_file, 'w') as o:
        writer = csv.writer(o)
        for seqid, taxid in raw_dict.items():
            gen_taxid = get_parent(taxid, taxonomy_dict, rank)
            writer.writerow([seqid, gen_taxid])



        

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--rank", default="genus")
    args = parser.parse_args()

    raw_dict = process_kraken('./data/kraken_results/output.txt')
    taxonomy_dict = load_NCBI(NCBI_taxonomy)
    process_raw(raw_dict, taxonomy_dict, args.rank, './data/kraken_results/output.csv')