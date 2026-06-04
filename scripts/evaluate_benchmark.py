import csv
import argparse


def load_abundance(file_tsv):
    """Read file TSV containing the abundances of each contig"""
    contig_ab = {}
    with open(file_tsv, 'r') as f:
        for line in f:
            columns = line.strip().split()
            if len(columns) >= 2:
                try:
                    contig_ab[columns[0]] = float(columns[1])
                except ValueError:
                    continue
    return contig_ab


def load_truth(truth_file):
    """Load the truth file containing the real abundance of each specie. Expected format (TSV): TaxID Relative_Abundance"""
    truth = {}
    with open(truth_file, 'r') as f:
        for line in f:
            columns = line.strip().split()
            if len(columns) >= 2:
                truth[columns[0]] = float(columns[1])

    # Normalize to percentages
    sum_truth = sum(truth.values())
    if sum_truth > 0:
        truth = {key: val / sum_truth for key, val in truth.items()}
    return truth


def map_bin_to_taxid(kraken_out, graphbin_out):
    """Assign TaxID to the GraphBin2 bins according to the Kraken seeds"""
    contig_to_taxid = {}
    with open(kraken_out, 'r') as f:
        for line in f:
            columns = line.strip().split()
            if len(columns) >= 3 and columns[0] == 'C':
                contig_to_taxid[columns[1]] = columns[2]

    bin_taxid_counts = {}
    with open(graphbin_out, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 2:
                contig_id = row[0]
                bin_id = row[1]
                if contig_id in contig_to_taxid:
                    taxid = contig_to_taxid[contig_id]
                    if bin_id not in bin_taxid_counts:
                        bin_taxid_counts[bin_id] = {}
                    bin_taxid_counts[bin_id][taxid] = bin_taxid_counts[bin_id].get(taxid, 0) + 1

    # Assign most frequent TaxID to the bin (majority vote)
    bin_to_taxid = {}
    for bin_id, counts in bin_taxid_counts.items():
        best_taxid = max(counts, key=counts.get)
        bin_to_taxid[bin_id] = best_taxid

    return bin_to_taxid

def compute_abundance_profile(assignment_file, contig_ab, format="kraken", map_bin=None):
    """Compute the total relative abundance for each TaxID"""
    # Sum abundances of single contigs to obtain total abundance of taxon/bin.
    profile = {}
    sample_total = 0.0

    with open(assignment_file, 'r') as f:
        if format == "kraken":
            for line in f:
                columns = line.strip().split()
                if len(columns) >= 3 and columns[0] == 'C':  # Only classified
                    contig_id = columns[1]
                    taxon_id = columns[2]

                    if contig_id in contig_ab:
                        val = contig_ab[contig_id]
                        profile[taxon_id] = profile.get(taxon_id, 0.0) + val
                        sample_total += val

        elif format == "graphbin":
            reader = csv.reader(f)
            for row in reader:
                if len(row) >= 2:
                    contig_id = row[0]
                    bin_id = row[1]

                    taxon_id = map_bin.get(bin_id, "Unknown") if map_bin else bin_id

                    if contig_id in contig_ab and taxon_id != "Unknown":
                        val = contig_ab[contig_id]
                        profile[taxon_id] = profile.get(taxon_id, 0.0) + val
                        sample_total += val

    # Normalize to percentages
    normalized_profile = {}
    if sample_total > 0:
        for key, val in profile.items():
            normalized_profile[key] = val / sample_total

    return normalized_profile

def compute_read_profile(file_path, format="kraken", map_bin=None):
    """Compute the relative abundance based on literal read counts"""
    profile = {}
    total_classified_reads = 0

    with open(file_path, 'r') as f:
        if format == "kraken":
            for line in f:
                columns = line.strip().split()
                if len(columns) >= 3 and columns[0] == 'C':
                    taxid = columns[2]
                    profile[taxid] = profile.get(taxid, 0) + 1
                    total_classified_reads += 1
        elif format == "graphbin":
            reader = csv.reader(f)
            next(reader, None)  # Skip header
            for row in reader:
                if len(row) >= 3:
                    bin_id = row[1]
                    taxid = map_bin.get(bin_id, "Unknown") if map_bin else bin_id
                    if taxid != "Unknown":
                        profile[taxid] = profile.get(taxid, 0) + 1
                        total_classified_reads += 1

    normalized_profile = {}
    if total_classified_reads > 0:
        for key, val in profile.items():
            normalized_profile[key] = val / total_classified_reads

    return normalized_profile, total_classified_reads

def compute_bray_curtis(estimation, real):
    """Compute Bray-Curtis distance (0 = Perfect, 1 = Completely wrong)"""
    common_taxa = set(estimation.keys()).union(set(real.keys()))
    sum_differences = 0.0
    total_sum = 0.0

    for taxon in common_taxa:
        estimated_val = estimation.get(taxon, 0.0)
        real_val = real.get(taxon, 0.0)
        sum_differences += abs(estimated_val - real_val)
        total_sum += (estimated_val + real_val)

    if total_sum == 0: return 1.0
    return sum_differences / total_sum

def analyze_results(kraken_out, graphbin_out, abundance_file, truth_file, kraken_reads, graphbin_reads):
    contig_ab = load_abundance(abundance_file)
    truth = load_truth(truth_file)

    print("=== ANALYSIS BASED ON CONTIGS ===")

    kraken_profile = compute_abundance_profile(kraken_out, contig_ab, format="kraken")
    bc_kraken = compute_bray_curtis(kraken_profile, truth)
    print("\nKRAKEN")
    print(f"- Identified taxa: {len(kraken_profile)}")
    print(f"- Bray-Curtis distance: {bc_kraken:.4f}")

    map_bin = None
    if graphbin_out:
        map_bin = map_bin_to_taxid(kraken_out, graphbin_out)
        graphbin_profile = compute_abundance_profile(graphbin_out, contig_ab, format="graphbin", map_bin=map_bin)
        bc_graphbin = compute_bray_curtis(graphbin_profile, truth)
        print("\nGRAPHBIN2")
        print(f"- Identified taxa after propagation: {len(graphbin_profile)}")
        print(f"- Bray-Curtis distance: {bc_graphbin:.4f}")

        print("\nCONCLUSION")
        if bc_kraken == 0.0 and bc_graphbin == 0.0:
            print("Both methods are perfect (distance 0.0).")
        elif bc_kraken == 0.0 and bc_graphbin > 0.0:
            print(
                f"Kraken is perfect (distance 0.0). GraphBin2 brings distance to {bc_graphbin:.4f}.")
        elif bc_graphbin < bc_kraken:
            improvement = (bc_kraken - bc_graphbin) / bc_kraken * 100
            print(f"GraphBin2 improves the prediction by {improvement:.1f}% compared to Kraken.")
        elif bc_graphbin > bc_kraken:
            worsening = (bc_graphbin - bc_kraken) / bc_kraken * 100
            print(f"GraphBin2 worsens the estimation by {worsening:.1f}% compared to Kraken.")
        else:
            print("GraphBin2 does not alterate the estimation compared to Kraken.")

    if kraken_reads and graphbin_reads and map_bin:
        print("\n=== ANALYSIS BASED ON READS ===")

        kr_profile, kr_total = compute_read_profile(kraken_reads, format="kraken")
        bc_kr = compute_bray_curtis(kr_profile, truth)
        print("\nKRAKEN")
        print(f"- Classified reads: {kr_total}")
        print(f"- Bray-Curtis distance: {bc_kr:.4f}")

        gb_profile, gb_total = compute_read_profile(graphbin_reads, format="graphbin", map_bin=map_bin)
        bc_gb = compute_bray_curtis(gb_profile, truth)
        print("\nGRAPHBIN")
        print(f"- Recovered reads: {gb_total}")
        print(f"- Bray-Curtis distance: {bc_gb:.4f}")

        diff_reads = gb_total - kr_total
        if diff_reads > 0:
            print(
                f"\nPerfect! GraphBin2 allows us to classify {diff_reads} reads more than Kraken.")
        else:
            print(
                f"\nAttention: Kraken classifies {abs(diff_reads)} reads more than the GraphBin2 method.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--kraken", required=True, help="Original output of Kraken on contigs")
    parser.add_argument("--abundance", required=True, help="File TSV with contigs' abundaces")
    parser.add_argument("--truth", required=True, help="File TSV with real abundances (TaxID vs Perc)")
    parser.add_argument("--graphbin", required=False, help="Output CSV of GraphBin2 (optional)")
    parser.add_argument("--kraken_reads", required=False, help="Original output of Kraken on raw reads")
    parser.add_argument("--graphbin_reads", required=False, help="Output CSV of reads_binning.py")

    args = parser.parse_args()
    analyze_results(args.kraken, args.graphbin, args.abundance, args.truth, args.kraken_reads, args.graphbin_reads)