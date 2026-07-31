import argparse
import csv
import gzip

def generate_binning(input_file, file_type, mapping_file, out_path, sample_id):

    with open(out_path, 'w') as out:
        # Write the .binning file following the CAMI format
        out.write("@Version:0.9.1\n")
        out.write(f"@SampleID:{sample_id}\n")
        out.write("@@SEQUENCEID\tBINID\tTAXID\n")

        if file_type == 'kraken':
            with open(input_file, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    # Take only classified reads
                    if len(parts) >= 3 and parts[0] == 'C':
                        seq_id = parts[1].split('/')[0]
                        tax_id = parts[2]
                        out.write(f"{seq_id}\t{tax_id}\t{tax_id}\t1\n")

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
                    seq_id = row['read_id']
                    bin_id = row['taxonomy']
                    # Ignore reads that GraphBin does not classify
                    if bin_id in bin_to_taxid and bin_id.lower() != "unclassified":
                        tax_id = bin_to_taxid[bin_id]
                        out.write(f"{seq_id}\t{tax_id}\t{tax_id}\t1\n")

        elif file_type == 'goldstandard':
            if input_file.endswith('.gz'):
                f = gzip.open(input_file, 'rt')
            else:
                f = open(input_file, 'rt')

            with f:
                for line in f:
                    # Skip headers of CAMI
                    if line.startswith('@') or line.startswith('#'):
                        continue
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        seq_id = parts[0]
                        bin_id = parts[1]
                        tax_id = parts[2]
                        out.write(f"{seq_id}\t{bin_id}\t{tax_id}\t1\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate CAMI .binning file from Kraken/GraphBin outputs")
    parser.add_argument("--input", required=True, help="Output reads file from Kraken/GraphBin")
    parser.add_argument("--type", choices=['kraken', 'graphbin', 'goldstandard'], required=True, help="Tool used")
    parser.add_argument("--mapping", required=False, help="CSV file for mapping bin -> taxid from processor.py")
    parser.add_argument("--output", required=True, help="Output .binning file")
    parser.add_argument("--sample-id", default="Sample_1", help="Sample ID for header CAMI")
    args = parser.parse_args()

    generate_binning(args.input, args.type, args.mapping, args.output, args.sample_id)