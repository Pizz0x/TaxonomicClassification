import argparse
import csv
from collections import defaultdict

# Standard CAMI ranks
RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']

def load_taxonomy(nodes_file, names_file):
    """Load NCBI taxonomic tree"""
    nodes, ranks, names = {}, {}, {}

    with open(nodes_file, 'r') as f:
        for line in f:
            parts = line.split('|')
            if len(parts) >= 3:
                taxid = parts[0].strip()
                parent = parts[1].strip()
                rank = parts[2].strip()
                nodes[taxid] = parent
                ranks[taxid] = rank

    with open(names_file, 'r') as f:
        for line in f:
            if 'scientific name' in line:
                parts = line.split('|')
                if len(parts) >= 2:
                    names[parts[0].strip()] = parts[1].strip()

    return nodes, ranks, names


def get_taxpath(taxid, nodes, ranks, names):
    """Build TAXPATH and TAXPATHSN required by CAMI format"""
    path_ids, path_names = [], []
    curr = taxid

    # Climb the tree until the root
    lineage = []
    while curr and curr != '1' and curr != '0' and curr in nodes:
        lineage.append(curr)
        if curr == nodes[curr]:
            break
        curr = nodes[curr]

    lineage.reverse()  # From root to leaf

    # Filter only standard ranks CAMI
    for t in lineage:
        if ranks.get(t) in RANKS or t == taxid:
            path_ids.append(t)
            path_names.append(names.get(t, f"Unknown_{t}"))

    return "|".join(path_ids), "|".join(path_names)


def generate_profile(input_file, file_type, mapping_file, nodes_file, names_file, out_path, sample_id, total_reads):

    if not total_reads or total_reads <= 0:
        print("Error: invalid total_reads value")
        return

    nodes, ranks, names = load_taxonomy(nodes_file, names_file)
    read_counts = defaultdict(int)

    # Load data basing on the tool format
    if file_type == 'kraken':
        with open(input_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    if parts[0] == 'C':
                        taxid = parts[2]
                        read_counts[taxid] += 1

    elif file_type == 'graphbin':
        bin_to_taxid = {}
        with open(mapping_file, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) >= 2:
                    bin_to_taxid[row[0]] = row[1]

        with open(input_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                bin_id = row['taxonomy']
                if bin_id in bin_to_taxid:
                    taxid = bin_to_taxid[bin_id]
                    read_counts[taxid] += 1

    # Propagate counts along the taxonomic tree
    cumulative_counts = defaultdict(int)
    for taxid, count in read_counts.items():
        curr = taxid
        while curr and curr != '1' and curr != '0' and curr in nodes:
            cumulative_counts[curr] += count
            if curr == nodes[curr]:
                break
            curr = nodes[curr]

    # Write the .profile file following the CAMI format
    with open(out_path, 'w') as out:
        out.write(f"@SampleID:{sample_id}\n")
        out.write("@Version:0.9.1\n")
        out.write("@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n")
        out.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")

        for taxid, count in cumulative_counts.items():
            rank = ranks.get(taxid, "unknown")
            if rank in RANKS or taxid in read_counts:
                percentage = (count / total_reads) * 100
                taxpath, taxpathsn = get_taxpath(taxid, nodes, ranks, names)
                out.write(f"{taxid}\t{rank}\t{taxpath}\t{taxpathsn}\t{percentage:.6f}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate CAMI .profile from Kraken or GraphBin")
    parser.add_argument("--input", required=True, help="Output reads file from the tool used")
    parser.add_argument("--type", required=True, choices=['kraken', 'graphbin'], help="Tool used")
    parser.add_argument("--mapping", required=False, help="CSV file for mapping bin -> taxid from processor.py")
    parser.add_argument("--nodes", required=True, help="Path to nodes.dmp from Kraken db")
    parser.add_argument("--names", required=True, help="Path to names.dmp from Kraken db")
    parser.add_argument("--output", required=True, help="Output .profile file")
    parser.add_argument("--sample-id", required=False, help="Sample ID for header CAMI")
    parser.add_argument("--total-reads", required=False, type=int, help="Total number of input reads")
    args = parser.parse_args()

    generate_profile(args.input, args.type, args.mapping, args.nodes, args.names, args.output, args.sample_id, args.total_reads)
