import shutil
import subprocess
import sys
from pathlib import Path

from Bio.PDB import PDBParser


def extract_sequence(pdb_file, model_id=0, chain_id=None):
    """
    Extract the sequence from a PDB file.

    Args:
        pdb_file (str): Path to the PDB file.
        model_id (int): Model ID to extract (default: 0, the first model).
                        Useful for NMR structures with multiple models.
        chain_id (str): Chain ID to extract (default: None, uses first chain only).
                        Set to 'all' to extract all chains.

    Returns:
        str: The extracted sequence as a string.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("RNA", pdb_file)
    sequence = []
    model = structure[model_id]

    if chain_id == "all":
        chains = model
    elif chain_id is not None:
        chains = [model[chain_id]]
    else:
        # Default: use first chain only
        chains = [next(iter(model))]

    for chain in chains:
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

    # Clean up invalid patterns like "()" - adjacent open/close with nothing between
    dot_bracket_str = fix_invalid_pairs(dot_bracket_str)

    print(dot_bracket_str)
    return dot_bracket_str if dot_bracket else None


def fix_invalid_pairs(dot_bracket):
    """
    Fix invalid dot-bracket patterns like "()" by converting them to "..".

    Args:
        dot_bracket (str): Dot-bracket notation string.

    Returns:
        str: Cleaned dot-bracket notation.
    """
    result = list(dot_bracket)
    changed = True

    # Repeatedly remove invalid "()" pairs until none remain
    while changed:
        changed = False
        i = 0
        while i < len(result) - 1:
            if result[i] == "(" and result[i + 1] == ")":
                result[i] = "."
                result[i + 1] = "."
                changed = True
            i += 1

    return "".join(result)


def write_fasta(output_path, header, sequence, dot_bracket):
    """
    Write the sequence and dot-bracket notation to a FASTA file.

    Args:
        output_path (str): Path to the output FASTA file.
        header (str): FASTA header (without >).
        sequence (str): RNA sequence.
        dot_bracket (str): Dot-bracket notation.
    """
    with open(output_path, "w") as fasta_file:
        fasta_file.write(f">{header}\n")
        fasta_file.write(f"{sequence}\n")
        fasta_file.write(f"{dot_bracket}\n")


def run_r2dt(fasta_file, output_name=None):
    """
    Run R2DT on the given FASTA file and copy the generated SVG file.

    Args:
        fasta_file (str): Path to the FASTA file.
        output_name (str): Optional custom name for output SVG (without extension).

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
        svg_destination_dir = input_folder / "svg"
        svg_destination_dir.mkdir(parents=True, exist_ok=True)
        if svg_source_dir.exists() and svg_source_dir.is_dir():
            for svg_file in svg_source_dir.glob("*.svg"):
                if output_name:
                    # Use custom output name
                    svg_destination = svg_destination_dir / f"{output_name}.svg"
                else:
                    svg_destination = svg_destination_dir / svg_file.name
                svg_file.replace(svg_destination)

        # Force delete the r2dt directory and its contents
        shutil.rmtree(r2dt_dir)

        return result.stdout
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Error: {e}")
        return None


def get_structure_info(pdb_file):
    """
    Get information about models and chains in a PDB file.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        tuple: (structure, list of (model_id, chain_id) tuples)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("RNA", pdb_file)

    model_chain_pairs = []
    for model in structure:
        for chain in model:
            # Check if chain has RNA residues
            has_rna = any(
                residue.id[0] == " " and residue.resname in ["A", "U", "G", "C"]
                for residue in chain
            )
            if has_rna:
                model_chain_pairs.append((model.id, chain.id))

    return structure, model_chain_pairs


def process_single(pdb_file, model_id=0, chain_id=None):
    """
    Process a single model/chain combination.

    Args:
        pdb_file (str): Path to the PDB file.
        model_id (int): Model ID to extract.
        chain_id (str): Chain ID to extract (None for first chain).

    Returns:
        bool: True if successful.
    """
    sequence = extract_sequence(pdb_file, model_id=model_id, chain_id=chain_id)
    if not sequence:
        print(f"  No sequence found for model {model_id}, chain {chain_id}")
        return False

    rnaview_output = run_rnaview(pdb_file)
    dot_bracket = parse_rnaview_output(rnaview_output, sequence)
    if not dot_bracket:
        print(f"  Error parsing RNAview output for model {model_id}, chain {chain_id}")
        return False

    # Generate output name
    pdb_path = Path(pdb_file)
    pdb_name = pdb_path.stem
    output_name = f"{pdb_name}_model{model_id}_chain{chain_id}"

    # Write FASTA
    fasta_path = pdb_path.parent / f"{output_name}.fasta"
    write_fasta(str(fasta_path), output_name, sequence, dot_bracket)

    print(f"  {output_name}: {len(sequence)} nt")

    # Run R2DT
    result = run_r2dt(str(fasta_path), output_name)

    # Clean up FASTA file
    if fasta_path.exists():
        fasta_path.unlink()

    return result is not None


def main():
    """
    Main function for processing a PDB file, generating RNA secondary structure,
    and running R2DT. Processes all models and chains.
    """
    if len(sys.argv) < 2:
        print("Usage: python3 rnaview.py <pdb_file> [--single]")
        print("  --single: Only process first model/chain (default: all)")
        sys.exit(1)

    pdb_file = sys.argv[1]
    single_mode = "--single" in sys.argv

    if single_mode:
        # Original behavior: first model, first chain
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

        pdb_path = Path(pdb_file)
        fasta_path = pdb_path.with_suffix(".fasta")
        write_fasta(str(fasta_path), pdb_path.stem, sequence, dot_bracket)
        r2dt_output = run_r2dt(str(fasta_path))
        print(r2dt_output)
    else:
        # Process all models and chains
        _, model_chain_pairs = get_structure_info(pdb_file)

        if not model_chain_pairs:
            print("No RNA chains found in PDB file.")
            sys.exit(1)

        # Check if multiple models (NMR) or just one
        unique_models = set(m for m, c in model_chain_pairs)
        unique_chains = set(c for m, c in model_chain_pairs if m == 0)

        print(
            f"Found {len(unique_models)} model(s) and {len(unique_chains)} RNA chain(s)"
        )

        # For NMR with many models, only process first model
        if len(unique_models) > 1:
            print(
                f"NMR structure detected ({len(unique_models)} models). Processing model 0 only."
            )
            model_chain_pairs = [(m, c) for m, c in model_chain_pairs if m == 0]

        success = 0
        for model_id, chain_id in model_chain_pairs:
            if process_single(pdb_file, model_id, chain_id):
                success += 1

        print(f"\nGenerated {success}/{len(model_chain_pairs)} images")


if __name__ == "__main__":
    main()
