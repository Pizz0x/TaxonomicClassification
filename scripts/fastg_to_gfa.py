
import argparse
import sys
import re


def fastg_to_gfa(fasta_file, gfa_file):
    """
    Convert FASTA/FASTG contigs to GFA format.
    Preserves contig names and extracts length information.
    Skips reverse complement entries (those with ' suffix).
    Extracts base node name from composite node names (e.g., NODE_1:NODE_1,NODE_2 -> NODE_1).
    """
    try:
        with open(fasta_file, 'r') as fi, open(gfa_file, 'w') as fo:
            name = None
            seq = []
            
            for line in fi:
                line = line.rstrip('\n')
                
                if line.startswith('>'):
                    # Write previous entry if exists and is not reverse complement
                    if name is not None and not name.endswith("'"):
                        sequence = ''.join(seq)
                        length = len(sequence)
                        fo.write(f"S\t{name}\t{sequence}\tLN:i:{length}\n")
                    
                    # Parse header (remove '>' and extract name)
                    header = line[1:]  # Remove '>'
                    # Extract just the contig name (first field before space or semicolon)
                    name = header.split()[0].rstrip(';')
                    # Extract base node name from composite node names (before any ':')
                    name = name.split(':')[0]
                    seq = []
                else:
                    # Accumulate sequence
                    seq.append(line.strip())
            
            # Write final entry if exists and is not reverse complement
            if name is not None and not name.endswith("'"):
                sequence = ''.join(seq)
                length = len(sequence)
                fo.write(f"S\t{name}\t{sequence}\tLN:i:{length}\n")
        
        print(f"Successfully converted {fasta_file} to {gfa_file}", file=sys.stderr)
        return 0
    
    except Exception as e:
        print(f"Error converting FASTA to GFA: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--input')
    parser.add_argument('--output')
    args = parser.parse_args()
    sys.exit(fastg_to_gfa(args.input, args.output))
