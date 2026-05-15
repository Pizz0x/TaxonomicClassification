import pysam, argparse

def assign_reads_with_filter(bam_input, csv_input, output_file):
    contig_to_tax = {}
    # create a dictionaty mapping each contigs to its taxonomy label
    with open(csv_input, "r") as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) == 2:
                contig_to_tax[parts[0]] = parts[1]

    # open the BAM file and the output CSV file, and write the header
    with pysam.AlignmentFile(bam_input, "rb") as bam, open(output_file, "w") as out:
        out.write("read_id,taxonomy,mapq\n")
        
        for read in bam:
            # check if read is mapped and has a mapping quality above the threshold
            if not read.is_unmapped:
                if not (read.is_secondary or read.is_supplementary):
                    if read.reference_name in contig_to_tax:
                        label = contig_to_tax[read.reference_name]
                        out.write(f"{read.query_name},{label},{read.mapping_quality}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment")
    parser.add_argument("--bin")
    parser.add_argument("--o")
    args = parser.parse_args()
    assign_reads_with_filter(args.alignment, args.bin, args.o)