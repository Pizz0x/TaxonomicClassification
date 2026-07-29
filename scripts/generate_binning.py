import argparse
import csv

def generate_binning(input_file, file_type, out_path, sample_id):

    with open(out_path, 'w') as out:
        # Write the .binning file following the CAMI format
        out.write("@Version:0.9.1\n")
        out.write(f"@SampleID:{sample_id}\n")
        out.write("@@SEQUENCEID\tBINID\tTAXID\n")

        if file_type == 'kraken':
            with open(input_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    # Take only classified reads
                    if len(parts) >= 3 and parts[0] == 'C':
                        seq_id = parts[1]
                        tax_id = parts[2]
                        out.write(f"{seq_id}\t{tax_id}\t{tax_id}\n")

        elif file_type == 'graphbin':
            with open(input_file, 'r') as f:
                reader = csv.reader(f)
                header = next(reader, None)  # Skip header
                for row in reader:
                    if len(row) >= 2:
                        seq_id = row[0].strip()
                        tax_id = row[1].strip()
                        # Ignore reads that GraphBin does not classify
                        if tax_id and tax_id.lower() != "unclassified":
                            out.write(f"{seq_id}\t{tax_id}\t{tax_id}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate CAMI .binning file from Kraken/GraphBin outputs")
    parser.add_argument("--input", required=True, help="Output reads file from Kraken/GraphBin")
    parser.add_argument("--type", choices=['kraken', 'graphbin'], required=True, help="Tool used")
    parser.add_argument("--output", required=True, help="Output .binning file")
    parser.add_argument("--sample-id", default="Sample_1", help="Sample ID for header CAMI")
    args = parser.parse_args()

    generate_binning(args.input, args.type, args.output, args.sample_id)