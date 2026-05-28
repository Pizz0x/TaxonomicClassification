<<<<<<< HEAD
import pysam, argparse, os, sys

"""
Assign reads to taxonomic bins using a contig-to-label CSV map.

Expected input formats:
- BAM: standard BAM file containing mapped reads
- CSV: two comma-separated columns: contig_id,label
  A header row is allowed and will be skipped if the first column contains "contig_id".

Output format:
- CSV with header read_id,taxonomy,mapq
"""

def file_exists(path, desc):
    if not path:
        raise FileNotFoundError(f"{desc} path is required")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{desc} not found: {path}")
    return path


def assign_reads_with_filter(bam_input, csv_input, output_file, min_mapq):
=======
import pysam, argparse

def assign_reads_with_filter(bam_input, csv_input, output_file):
>>>>>>> df16b5496842f1b4ed9efec35bddacaf51103357
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
<<<<<<< HEAD
            if not read.is_unmapped and read.mapping_quality >= min_mapq:
=======
            if not read.is_unmapped:
>>>>>>> df16b5496842f1b4ed9efec35bddacaf51103357
                if not (read.is_secondary or read.is_supplementary):
                    if read.reference_name in contig_to_tax:
                        label = contig_to_tax[read.reference_name]
                        out.write(f"{read.query_name},{label},{read.mapping_quality}\n")

if __name__ == "__main__":
<<<<<<< HEAD
    parser = argparse.ArgumentParser(description="Assign mapped reads to bins using a contig label map.")
    parser.add_argument("--alignment", required=True,
                        help="Input BAM file containing aligned reads.")
    parser.add_argument("--bin", required=True,
                        help="CSV file mapping contig IDs to bin labels. Format: contig_id,label")
    parser.add_argument("--o", required=True,
                        help="Output CSV file containing read_id,taxonomy,mapq")
    parser.add_argument("--min-mapq", type=int, default=0,
                        help="Minimum mapping quality to include a read.")
    args = parser.parse_args()
    assign_reads_with_filter(args.alignment, args.bin, args.o, args.min_mapq)



=======
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment")
    parser.add_argument("--bin")
    parser.add_argument("--o")
    args = parser.parse_args()
    assign_reads_with_filter(args.alignment, args.bin, args.o)
>>>>>>> df16b5496842f1b4ed9efec35bddacaf51103357
