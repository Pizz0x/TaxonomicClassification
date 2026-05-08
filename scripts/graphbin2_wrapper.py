#!/usr/bin/env python3
"""
Wrapper script for GraphBin2 that handles name mapping between k141_* contig IDs
(from kraken output and contigs file) and NODE_* IDs (from GFA file).
"""

import os
import sys
import subprocess
import csv
import tempfile
from collections import OrderedDict


def create_k141_to_node_mapping(contigs_file, gfa_file):
    """
    Create a mapping from k141_* contig IDs to NODE_* IDs based on sequence content.
    """
    # Read k141 sequences
    k141_seqs = {}
    with open(contigs_file, 'r') as f:
        current_name = None
        current_seq = []
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if current_name is not None:
                    seq = ''.join(current_seq)
                    k141_seqs[current_name] = seq
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_name is not None:
            seq = ''.join(current_seq)
            k141_seqs[current_name] = seq
    
    # Read NODE sequences from GFA
    node_seqs = {}
    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('S'):
                parts = line.rstrip('\n').split('\t')
                node_name = parts[1]
                sequence = parts[2]
                # Only add forward nodes (not reverse complements)
                if not node_name.endswith("'"):
                    node_seqs[node_name] = sequence
    
    # Create mapping based on sequence content
    mapping = {}
    seq_to_node = {seq: node for node, seq in node_seqs.items()}
    
    for k141_id, k141_seq in k141_seqs.items():
        if k141_seq in seq_to_node:
            node_id = seq_to_node[k141_seq]
            mapping[k141_id] = node_id
            print(f"Mapping (by sequence): {k141_id} -> {node_id}", file=sys.stderr)
        else:
            print(f"Warning: Could not find match for {k141_id} in GFA", file=sys.stderr)
    
    return mapping


def create_node_contigs_file(contigs_file, mapping, gfa_file):
    """
    Create a version of the contigs file using NODE names instead of k141 names.
    """
    # Create temp file with NODE names
    temp_fd, temp_path = tempfile.mkstemp(suffix='.fa', text=True)
    with os.fdopen(temp_fd, 'w') as f:
        # Read k141 contigs in order and write with NODE names
        with open(contigs_file, 'r') as infile:
            current_name = None
            current_seq = []
            current_header = None
            for line in infile:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if current_name is not None and current_name in mapping:
                        seq = ''.join(current_seq)
                        # Use NODE name but preserve other header fields
                        node_name = mapping[current_name]
                        if current_header:
                            header_parts = current_header[1:].split()
                            # Replace contig ID with NODE name
                            header_parts[0] = node_name
                            new_header = ' '.join(header_parts)
                        else:
                            new_header = node_name
                        f.write(f">{new_header}\n{seq}\n")
                    
                    current_header = line
                    header_fields = line[1:].split()
                    current_name = header_fields[0]
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            
            if current_name is not None and current_name in mapping:
                seq = ''.join(current_seq)
                node_name = mapping[current_name]
                if current_header:
                    header_parts = current_header[1:].split()
                    header_parts[0] = node_name
                    new_header = ' '.join(header_parts)
                else:
                    new_header = node_name
                f.write(f">{new_header}\n{seq}\n")
    
    return temp_path


def run_graphbin2(args):
    """
    Run GraphBin2 with NODE-based naming and convert output to use k141_* names.
    """
    # Parse graphbin2 command-line arguments
    contigs_file = None
    gfa_file = None
    binned_file = None
    output_path = None
    remaining_args = []
    
    i = 0
    while i < len(args):
        if args[i] == '--contigs' and i + 1 < len(args):
            contigs_file = args[i + 1]
            i += 2
        elif args[i] == '--graph' and i + 1 < len(args):
            gfa_file = args[i + 1]
            i += 2
        elif args[i] == '--binned' and i + 1 < len(args):
            binned_file = args[i + 1]
            i += 2
        elif args[i] == '--output' and i + 1 < len(args):
            output_path = args[i + 1]
            i += 2
        else:
            remaining_args.append(args[i])
            i += 1
    
    if not (contigs_file and gfa_file and binned_file and output_path):
        print("Error: Missing required arguments", file=sys.stderr)
        sys.exit(1)
    
    # Create mapping
    mapping = create_k141_to_node_mapping(contigs_file, gfa_file)
    
    # Create temp files for NODE-named versions
    temp_node_contigs = create_node_contigs_file(contigs_file, mapping, gfa_file)
    
    temp_binned_file = binned_file + '.temp_node'
    with open(binned_file, 'r') as infile, open(temp_binned_file, 'w') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        for row in reader:
            if row[0] in mapping:
                row[0] = mapping[row[0]]
            writer.writerow(row)
    
    # Build GraphBin2 command with all arguments
    graphbin2_args = ['graphbin2'] + remaining_args + [
        '--graph', gfa_file,
        '--contigs', temp_node_contigs,
        '--binned', temp_binned_file,
        '--output', output_path
    ]
    
    print(f"Running: {' '.join(graphbin2_args)}", file=sys.stderr)
    result = subprocess.run(graphbin2_args)
    
    if result.returncode != 0:
        print("GraphBin2 failed", file=sys.stderr)
        # Clean up
        if os.path.exists(temp_node_contigs):
            os.remove(temp_node_contigs)
        if os.path.exists(temp_binned_file):
            os.remove(temp_binned_file)
        sys.exit(result.returncode)
    
    # Convert output back to k141_* names
    output_csv = os.path.join(output_path, 'graphbin2_output.csv')
    if os.path.exists(output_csv):
        temp_output = output_csv + '.temp'
        reverse_mapping = {v: k for k, v in mapping.items()}
        
        with open(output_csv, 'r') as infile, open(temp_output, 'w') as outfile:
            reader = csv.reader(infile)
            writer = csv.writer(outfile)
            for row in reader:
                if row[0] in reverse_mapping:
                    row[0] = reverse_mapping[row[0]]
                writer.writerow(row)
        
        os.replace(temp_output, output_csv)
        print(f"Converted output back to k141_* names", file=sys.stderr)
    
    # Clean up temp files
    if os.path.exists(temp_node_contigs):
        os.remove(temp_node_contigs)
    if os.path.exists(temp_binned_file):
        os.remove(temp_binned_file)
    
    return result.returncode


if __name__ == '__main__':
    sys.exit(run_graphbin2(sys.argv[1:]))
