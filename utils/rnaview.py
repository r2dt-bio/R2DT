import sys
import subprocess
from Bio import SeqIO
from Bio.PDB import PDBParser
import re
import os
from pathlib import Path
import shutil

def extract_sequence(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('RNA', pdb_file)
    sequence = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':
                    sequence.append(residue.resname)
    print("".join(sequence))
    return ''.join(sequence)


def run_rnaview(pdb_file):
    try:
        result = subprocess.run(['rnaview', pdb_file], capture_output=True, text=True)
        result.check_returncode()
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running RNAview: {e}")
        sys.exit(1)

def parse_rnaview_output(output, sequence):
    canonical_pairs = {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}
    dot_bracket = ['.'] * len(sequence)
    pairings = {}
    lines = output.splitlines()
    
    for line in lines:
        if line.startswith('pair'):
            parts = line.split()
            if len(parts) >= 5:
                pos1 = int(parts[2]) - 1  # Convert to 0-based index
                pos2 = int(parts[3]) - 1  # Convert to 0-based index
                if pos1 < len(sequence) and pos2 < len(sequence):
                    pair = sequence[pos1] + sequence[pos2]
                    if pair in canonical_pairs:
                        pairings[pos1] = pos2
                        pairings[pos2] = pos1
    
    for pos1, pos2 in pairings.items():
        if pos1 < pos2:
            dot_bracket[pos1] = '('
            dot_bracket[pos2] = ')'
    print("".join(dot_bracket))
    return ''.join(dot_bracket) if dot_bracket else None

def _extract_sequence(pdb_file):
    with open(pdb_file, 'r') as file:
        for record in SeqIO.parse(file, 'pdb-atom'):
            return str(record.seq)
    return None

def write_fasta(pdb_file, sequence, dot_bracket):
    fasta_filename = pdb_file.replace(".pdb", ".fasta")
    base_filename = os.path.basename(pdb_file).replace(".pdb", "")
    with open(fasta_filename, 'w') as fasta_file:
        fasta_file.write(f">{base_filename}\n")
        fasta_file.write(f"{sequence}\n")
        fasta_file.write(f"{dot_bracket}\n")

def run_r2dt(fasta_file):
    fasta_path = Path(fasta_file)
    input_folder = fasta_path.parent
    filename_stem = fasta_path.stem
    r2dt_dir = input_folder / 'r2dt'
    
    # Check if the directory exists, if not, create it
    if not r2dt_dir.exists():
        r2dt_dir.mkdir(parents=True)
    
    try:
        output_path = r2dt_dir / filename_stem
        result = subprocess.run(['r2dt.py', 'templatefree', str(fasta_path), str(output_path)], capture_output=True, text=True)
        result.check_returncode()
        # Copy the *.svg file before deletion
        svg_source_dir = r2dt_dir / filename_stem / 'results' / 'svg'
        svg_destination_dir = input_folder
        if svg_source_dir.exists() and svg_source_dir.is_dir():
            for svg_file in svg_source_dir.glob('*.svg'):
                svg_destination = svg_destination_dir / svg_file.name
                svg_file.replace(svg_destination)
        
        # Force delete the r2dt directory and its contents
        shutil.rmtree(r2dt_dir)
        
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running R2DT: {e}")
        sys.exit(1)
    except FileNotFoundError as e:
        print(f"File not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 rnaview.py <pdb_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    sequence = extract_sequence(pdb_file)
    if not sequence:
        print("Error extracting sequence from PDB file.")
        sys.exit(1)

    rnaview_output = run_rnaview(pdb_file)
    dot_bracket = parse_rnaview_output(rnaview_output, sequence)
    if not dot_bracket:
        print("Error parsing RNAview output.")
        sys.exit(1)

    print(f"Filename: {pdb_file}")
    print(f"PDB Sequence: {sequence}")
    print(f"Secondary Structure (Dot-Bracket Notation): {dot_bracket}")

    write_fasta(pdb_file, sequence, dot_bracket)
    fasta_file = pdb_file.replace(".pdb", ".fasta")
    r2dt_output = run_r2dt(fasta_file)
    print(r2dt_output)


if __name__ == "__main__":
    main()