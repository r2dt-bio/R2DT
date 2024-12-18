import os
import shutil
import subprocess
import sys
from pathlib import Path

from Bio.PDB import PDBParser


def extract_sequence(pdb_file):
    """
    Extract the sequence from a PDB file.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        str: The extracted sequence as a string.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("RNA", pdb_file)
    sequence = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":
                    sequence.append(residue.resname)
    print("".join(sequence))
    return "".join(sequence)


def run_rnaview(pdb_file):
    """
    Run RNAview on a given PDB file.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        str: RNAview output.
    """
    try:
        result = subprocess.run(
            ["rnaview", pdb_file], capture_output=True, text=True, check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running RNAview: {e}")
        sys.exit(1)


def parse_rnaview_output(output, sequence):
    """
    Parse RNAview output and convert it to dot-bracket notation.

    Args:
        output (str): RNAview output.
        sequence (str): RNA sequence.

    Returns:
        str: Dot-bracket notation or None if parsing fails.
    """
    canonical_pairs = {"AU", "UA", "GC", "CG", "GU", "UG"}
    dot_bracket = ["."] * len(sequence)
    pairings = {}  # Dictionary to hold sets of pairings
    lines = output.splitlines()

    # Step 1: Collect all pairings
    for line in lines:
        if not line.startswith("pair"):
            continue  # Skip lines that do not start with "pair"

        parts = line.split()
        if len(parts) >= 5:
            pos1 = int(parts[2]) - 1  # Convert to 0-based index
            pos2 = int(parts[3]) - 1  # Convert to 0-based index

            if pos1 < len(sequence) and pos2 < len(sequence):
                pair = sequence[pos1] + sequence[pos2]
                if pair in canonical_pairs:
                    # Initialize sets if not already
                    if pos1 not in pairings:
                        pairings[pos1] = set()
                    if pos2 not in pairings:
                        pairings[pos2] = set()

                    # Add the pairings
                    pairings[pos1].add(pos2)
                    pairings[pos2].add(pos1)

    # Step 2: Identify bases that pair with more than one base
    multi_paired_bases = {
        pos for pos, partners in pairings.items() if len(partners) > 1
    }

    # Step 3: Construct the dot-bracket notation
    for pos1, partners in pairings.items():
        # Skip if base is multi-paired
        if pos1 in multi_paired_bases:
            continue
        # Get the single partner
        pos2 = next(iter(partners))
        # Skip if partner is multi-paired or already assigned
        if (
            pos2 in multi_paired_bases
            or dot_bracket[pos1] != "."
            or dot_bracket[pos2] != "."
        ):
            continue
        if pos1 < pos2:
            dot_bracket[pos1] = "("
            dot_bracket[pos2] = ")"
        else:
            dot_bracket[pos2] = "("
            dot_bracket[pos1] = ")"

    # Output the dot-bracket notation
    dot_bracket_str = "".join(dot_bracket)
    print(dot_bracket_str)
    return dot_bracket_str if dot_bracket else None


def write_fasta(pdb_file, sequence, dot_bracket):
    """
    Write the sequence and dot-bracket notation to a FASTA file.

    Args:
        pdb_file (str): Path to the PDB file.
        sequence (str): RNA sequence.
        dot_bracket (str): Dot-bracket notation.
    """
    fasta_filename = pdb_file.replace(".pdb", ".fasta")
    base_filename = os.path.basename(pdb_file).replace(".pdb", "")
    with open(fasta_filename, "w") as fasta_file:
        fasta_file.write(f">{base_filename}\n")
        fasta_file.write(f"{sequence}\n")
        fasta_file.write(f"{dot_bracket}\n")


def run_r2dt(fasta_file):
    """
    Run R2DT on the given FASTA file and copy the generated SVG file.

    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        str: R2DT output.
    """
    fasta_path = Path(fasta_file)
    input_folder = fasta_path.parent
    filename_stem = fasta_path.stem
    r2dt_dir = input_folder / "r2dt"

    # Check if the directory exists, if not, create it
    if not r2dt_dir.exists():
        r2dt_dir.mkdir(parents=True)

    try:
        output_path = r2dt_dir / filename_stem
        result = subprocess.run(
            ["r2dt.py", "templatefree", str(fasta_path), str(output_path), "--quiet"],
            capture_output=True,
            text=True,
            check=True,
        )
        # Copy the *.svg file before deletion
        svg_source_dir = r2dt_dir / filename_stem / "results" / "svg"
        svg_destination_dir = input_folder
        if svg_source_dir.exists() and svg_source_dir.is_dir():
            for svg_file in svg_source_dir.glob("*.svg"):
                svg_destination = svg_destination_dir / svg_file.name
                svg_file.replace(svg_destination)

        # Force delete the r2dt directory and its contents
        shutil.rmtree(r2dt_dir)

        return result.stdout
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Error: {e}")
        sys.exit(1)


def main():
    """
    Main function for processing a PDB file, generating RNA secondary structure,
    and running R2DT.
    """
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
