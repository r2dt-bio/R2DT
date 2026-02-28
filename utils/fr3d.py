# pylint: disable=too-many-lines
"""
FR3D-python wrapper for extracting RNA secondary structure from 3D structures.

This module provides functions to run FR3D-python's NA_pairwise_interactions
to extract basepair annotations, and convert them to dot-bracket notation.

Copyright [2009-present] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from utils.pdb_fetch import DecompressedStructureFile, open_structure_file


def get_structure_info(structure_file: str) -> Dict:
    """
    Get information about chains and models in a structure file.

    Args:
        structure_file: Path to the structure file (.cif or .pdb).

    Returns:
        Dict with 'models' (list of model IDs) and 'chains'
        (dict of model -> list of chain IDs with RNA).
    """
    structure_path = Path(structure_file)
    info = {"models": [], "chains": {}, "format": structure_path.suffix.lower()}

    if structure_path.suffix.lower() == ".cif":
        info = _get_cif_structure_info(str(structure_path))
    else:
        info = _get_pdb_structure_info(str(structure_path))

    return info


def _get_cif_structure_info(cif_file: str) -> Dict:
    """Get structure info from mmCIF file using FR3D."""
    info = {"models": [], "chains": {}, "format": ".cif"}
    try:
        from fr3d.cif.reader import Cif  # pylint: disable=import-outside-toplevel

        with open_structure_file(Path(cif_file), "r") as f:
            structure = Cif(f).structure()

        # Get all RNA nucleotides
        bases = list(structure.residues(type=["RNA linking"]))
        if not bases:
            bases = list(structure.residues(type=["DNA linking"]))

        if not bases:
            return info

        # Collect unique models and chains
        models_chains = defaultdict(set)
        for base in bases:
            models_chains[base.model].add(base.chain)

        info["models"] = sorted(models_chains.keys())
        info["chains"] = {m: sorted(list(c)) for m, c in models_chains.items()}

    except ImportError:
        print("FR3D not installed.")
    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error reading CIF structure info: {e}")

    return info


# pylint: disable=too-many-nested-blocks
def _get_pdb_structure_info(pdb_file: str) -> Dict:
    """Get structure info from PDB file using BioPython."""
    info = {"models": [], "chains": {}, "format": ".pdb"}
    try:
        from Bio.PDB import PDBParser  # pylint: disable=import-outside-toplevel

        parser = PDBParser(QUIET=True)
        with DecompressedStructureFile(Path(pdb_file)) as decompressed_path:
            structure = parser.get_structure("RNA", str(decompressed_path))

        for model in structure:
            model_id = model.id
            chains_with_rna = []

            for chain in model:
                has_rna = False
                for residue in chain:
                    if residue.id[0] == " ":
                        resname = residue.resname.strip()
                        if resname in ["A", "C", "G", "U", "DA", "DC", "DG", "DT"]:
                            has_rna = True
                            break
                        if _get_parent_base(resname):
                            has_rna = True
                            break
                if has_rna:
                    chains_with_rna.append(chain.id)

            if chains_with_rna:
                info["models"].append(model_id)
                info["chains"][model_id] = sorted(chains_with_rna)

    except ImportError:
        print("BioPython not installed.")
    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error reading PDB structure info: {e}")

    return info


def run_fr3d(
    structure_file: str,
    output_dir: str,
    category: str = "basepair",
    quiet: bool = False,
) -> Optional[str]:
    """
    Run FR3D NA_pairwise_interactions.py on a structure file.

    FR3D can process both .cif and .pdb files directly.
    Supports gzip-compressed files (.gz).

    Args:
        structure_file: Path to the structure file (.cif or .pdb, optionally .gz).
        output_dir: Directory to write output files.
        category: Interaction category to annotate (default: "basepair").
        quiet: If True, suppress FR3D progress output.

    Returns:
        Path to the basepairs output file, or None if failed.
    """
    structure_path = Path(structure_file)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Extract PDB ID from filename (handle .gz extension)
    stem = structure_path.stem
    if stem.endswith(".cif") or stem.endswith(".pdb"):
        pdb_id = stem.rsplit(".", 1)[0]
    else:
        pdb_id = stem

    # Use DecompressedStructureFile to handle .gz files
    with DecompressedStructureFile(structure_path) as decompressed_path:
        try:
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "fr3d.classifiers.NA_pairwise_interactions",
                    "-i",
                    str(decompressed_path.parent),
                    "-o",
                    str(output_path),
                    "-c",
                    category,
                    decompressed_path.stem,
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            if not quiet:
                if result.stdout:
                    print(result.stdout)
                if result.stderr:
                    print(result.stderr, file=sys.stderr)

            # FR3D creates output file named {stem}_basepair.txt
            output_file = output_path / f"{decompressed_path.stem}_basepair.txt"
            if output_file.exists():
                # Rename to use original pdb_id if different
                if decompressed_path.stem != pdb_id:
                    final_output = output_path / f"{pdb_id}_basepair.txt"
                    output_file.rename(final_output)
                    return str(final_output)
                return str(output_file)
            return None

        except subprocess.CalledProcessError as e:
            if not quiet:
                print(f"Error running FR3D: {e}")
                print(f"stdout: {e.stdout}")
                print(f"stderr: {e.stderr}")
            return None


# pylint: disable=too-many-branches
def extract_sequence_from_cif(
    cif_file: str,
    chain_id: Optional[str] = None,
    model_id: Optional[int] = None,
    quiet: bool = False,
) -> Tuple[str, Dict[str, int]]:
    """
    Extract RNA sequence from a mmCIF file using FR3D's structure loading.

    Supports gzip-compressed files (.gz).

    Args:
        cif_file: Path to the mmCIF file.
        chain_id: Optional chain ID to extract. If None, extracts first RNA chain.
        model_id: Optional model ID to extract. If None, extracts first model only.
        quiet: If True, suppress error messages.

    Returns:
        Tuple of (sequence string, mapping of unit_id to sequence position).
    """
    try:
        from fr3d.cif.reader import Cif  # pylint: disable=import-outside-toplevel

        with open_structure_file(Path(cif_file), "r") as f:
            structure = Cif(f).structure()

        # Get RNA nucleotides
        bases = list(structure.residues(type=["RNA linking"]))

        if not bases:
            # Try DNA if no RNA found
            bases = list(structure.residues(type=["DNA linking"]))

        if not bases:
            return "", {}

        # Filter by model - use first model if not specified
        if model_id is not None:
            bases = [b for b in bases if b.model == model_id]
        else:
            # Use first model found
            first_model = bases[0].model if bases else None
            if first_model is not None:
                bases = [b for b in bases if b.model == first_model]

        # Filter by chain if specified
        if chain_id:
            bases = [b for b in bases if b.chain == chain_id]
        else:
            # Use first chain found (after model filtering)
            first_chain = bases[0].chain if bases else None
            if first_chain:
                bases = [b for b in bases if b.chain == first_chain]

        # Sort by index
        bases = sorted(bases, key=lambda b: (b.model, b.chain, b.index))

        # Build sequence and mapping
        sequence = []
        unit_id_to_position = {}

        for i, base in enumerate(bases):
            # Map modified nucleotides to standard bases
            seq = base.sequence
            if seq in ["A", "C", "G", "U"]:
                sequence.append(seq)
            elif seq in ["DA", "DC", "DG", "DT"]:
                sequence.append(seq[1] if seq != "DT" else "T")
            else:
                # Try to map modified nucleotide to parent
                parent = _get_parent_base(seq)
                sequence.append(parent if parent else "N")

            unit_id_to_position[base.unit_id()] = i

        return "".join(sequence), unit_id_to_position

    except ImportError:
        if not quiet:
            print("FR3D not installed. Please install fr3d-python.")
        return "", {}
    except Exception as e:  # pylint: disable=broad-exception-caught
        if not quiet:
            print(f"Error extracting sequence from CIF: {e}")
        return "", {}


# pylint: disable=too-many-branches,too-many-statements
def extract_sequence_from_pdb(
    pdb_file: str, chain_id: Optional[str] = None, quiet: bool = False
) -> Tuple[str, Dict[str, int]]:
    """
    Extract RNA sequence from a PDB file using BioPython.

    Supports gzip-compressed files (.gz).

    Args:
        pdb_file: Path to the PDB file.
        chain_id: Optional chain ID to extract. If None, extracts first RNA chain.
        quiet: If True, suppress error messages.

    Returns:
        Tuple of (sequence string, mapping of residue_key to sequence position).
    """
    try:
        from Bio.PDB import PDBParser  # pylint: disable=import-outside-toplevel

        parser = PDBParser(QUIET=True)
        with DecompressedStructureFile(Path(pdb_file)) as decompressed_path:
            structure = parser.get_structure("RNA", str(decompressed_path))

            # Get first model
            model = structure[0]

            # Collect RNA residues
            # When no chain is specified, prefer true RNA chains (A/C/G/U)
            # over DNA chains (DA/DC/DG/DT)
            rna_residues = []
            selected_chain = chain_id

            if not selected_chain:
                # Two-pass: first look for RNA, then fall back to DNA
                for prefer_rna in (True, False):
                    for chain in model:
                        has_target = False
                        for residue in chain:
                            if residue.id[0] == " ":
                                resname = residue.resname.strip()
                                if prefer_rna and resname in (
                                    "A",
                                    "C",
                                    "G",
                                    "U",
                                ):
                                    has_target = True
                                    break
                                if prefer_rna and _get_parent_base(resname):
                                    has_target = True
                                    break
                                if not prefer_rna and resname in (
                                    "DA",
                                    "DC",
                                    "DG",
                                    "DT",
                                ):
                                    has_target = True
                                    break
                        if has_target:
                            selected_chain = chain.id
                            break
                    if selected_chain:
                        break

            if not selected_chain:
                return "", {}

            for chain in model:
                if chain.id != selected_chain:
                    continue
                for residue in chain:
                    if residue.id[0] == " ":  # Standard residue
                        resname = residue.resname.strip()
                        if resname in ["A", "C", "G", "U", "DA", "DC", "DG", "DT"]:
                            rna_residues.append((chain.id, residue))
                        elif _get_parent_base(resname):
                            rna_residues.append((chain.id, residue))

            if not rna_residues:
                return "", {}

            # Sort by chain and residue number
            rna_residues.sort(key=lambda x: (x[0], x[1].id[1]))

            # Build sequence and mapping
            sequence = []
            residue_to_position = {}

            # Get pdb_id from original file, handling .gz extension
            pdb_path = Path(pdb_file)
            pdb_id = pdb_path.stem
            if pdb_id.endswith(".pdb") or pdb_id.endswith(".cif"):
                pdb_id = pdb_id.rsplit(".", 1)[0]

            for i, (chain, residue) in enumerate(rna_residues):
                resname = residue.resname.strip()
                if resname in ["A", "C", "G", "U"]:
                    sequence.append(resname)
                elif resname in ["DA", "DC", "DG", "DT"]:
                    sequence.append(resname[1] if resname != "DT" else "T")
                else:
                    parent = _get_parent_base(resname)
                    sequence.append(parent if parent else "N")

                # Create a key similar to FR3D unit_id format
                res_num = residue.id[1]
                key = f"{pdb_id}|1|{chain}|{resname}|{res_num}"
                residue_to_position[key] = i

            return "".join(sequence), residue_to_position

    except ImportError:
        if not quiet:
            print("BioPython not installed.")
        return "", {}
    except Exception as e:  # pylint: disable=broad-exception-caught
        if not quiet:
            print(f"Error extracting sequence from PDB: {e}")
        return "", {}


def _get_parent_base(modified_nt: str) -> Optional[str]:
    """
    Map modified nucleotide to parent base.

    Args:
        modified_nt: Modified nucleotide code.

    Returns:
        Parent base (A, C, G, or U) or None if unknown.
    """
    # Common modified nucleotides mapping
    modifications = {
        # Adenine derivatives
        "1MA": "A",
        "2MA": "A",
        "6MA": "A",
        "MIA": "A",
        "T6A": "A",
        # Cytosine derivatives
        "5MC": "C",
        "OMC": "C",
        "S4C": "C",
        # Guanine derivatives
        "1MG": "G",
        "2MG": "G",
        "7MG": "G",
        "M2G": "G",
        "OMG": "G",
        "YG": "G",
        # Uracil derivatives
        "5MU": "U",
        "2MU": "U",
        "4SU": "U",
        "H2U": "U",
        "PSU": "U",
        "OMU": "U",
        "5BU": "U",
        "UR3": "U",
        # Pseudouridine
        "P": "U",
    }
    return modifications.get(modified_nt.upper())


# ============================================================================
# Full-sequence extraction (SEQRES / entity_poly) for missing-residue support
# ============================================================================

# Standard three-letter RNA codes → one-letter
_THREE_TO_ONE = {
    "A": "A",
    "C": "C",
    "G": "G",
    "U": "U",
    "I": "I",
    "DA": "A",
    "DC": "C",
    "DG": "G",
    "DT": "T",
}


def _resname_to_one_letter(resname: str) -> str:
    """Convert a three-letter (or one-letter) residue name to single char."""
    resname = resname.strip().upper()
    if resname in _THREE_TO_ONE:
        return _THREE_TO_ONE[resname]
    parent = _get_parent_base(resname)
    return parent if parent else "N"


def get_full_sequence_from_pdb(  # pylint: disable=too-many-statements
    pdb_file: str, chain_id: Optional[str] = None
) -> Tuple[str, List[bool], str]:
    """
    Extract the full deposited sequence and resolved-residue mask from a PDB file.

    Parses SEQRES records for the complete sequence and compares with ATOM
    records to identify unresolved (missing) residue positions.

    Args:
        pdb_file: Path to PDB file (may be gzip-compressed).
        chain_id: Chain ID to extract.  If ``None`` the first RNA chain
            (based on ATOM records) is used and then the corresponding
            SEQRES chain is looked up.

    Returns:
        Tuple of ``(full_sequence, resolved_mask, used_chain_id)`` where

        * *full_sequence* is the deposited sequence from SEQRES,
        * *resolved_mask* is a list of booleans (``True`` = resolved),
        * *used_chain_id* is the chain that was actually used.

        Returns ``("", [], "")`` on failure.
    """
    try:
        from Bio.PDB import PDBParser  # pylint: disable=import-outside-toplevel

        parser = PDBParser(QUIET=True)

        with DecompressedStructureFile(Path(pdb_file)) as decompressed_path:
            # --- 1. Parse SEQRES lines manually (BioPython gives 'X' for RNA)
            seqres_by_chain: Dict[str, List[str]] = defaultdict(list)
            with open(str(decompressed_path), "r") as fh:
                for line in fh:
                    if line.startswith("SEQRES"):
                        cid = line[11].strip()
                        residues = line[19:].split()
                        seqres_by_chain[cid].extend(residues)

            # --- 2. Parse ATOM records for resolved residue numbers
            structure = parser.get_structure("RNA", str(decompressed_path))
            model = structure[0]

            # Determine chain if not provided
            # Prefer RNA chains (A/C/G/U) over DNA (DA/DC/DG/DT)
            if chain_id is None:
                dna_chain = None
                for chain_obj in model:
                    for residue in chain_obj:
                        if residue.id[0] == " ":
                            resname = residue.resname.strip()
                            if resname in ("A", "C", "G", "U") or _get_parent_base(
                                resname
                            ):
                                chain_id = chain_obj.id
                                break
                            if (
                                resname
                                in (
                                    "DA",
                                    "DC",
                                    "DG",
                                    "DT",
                                )
                                and dna_chain is None
                            ):
                                dna_chain = chain_obj.id
                    if chain_id:
                        break
                if chain_id is None:
                    chain_id = dna_chain

            if not chain_id:
                return "", [], ""

            # Collect resolved residue numbers for the chain
            resolved_resnums: Set[int] = set()
            if chain_id in [c.id for c in model]:
                for residue in model[chain_id]:
                    if residue.id[0] == " ":
                        resname = residue.resname.strip()
                        if resname in (
                            "A",
                            "C",
                            "G",
                            "U",
                            "DA",
                            "DC",
                            "DG",
                            "DT",
                        ) or _get_parent_base(resname):
                            resolved_resnums.add(residue.id[1])

            # --- 2b. Parse REMARK 465 for missing residue numbers
            missing_resnums: Set[int] = set()
            with open(str(decompressed_path), "r") as fh:
                for line in fh:
                    if line.startswith("REMARK 465") and len(line) > 20:
                        parts = line.split()
                        # Format: REMARK 465   resname chain resnum
                        if len(parts) >= 5:
                            try:
                                r465_chain = parts[3]
                                r465_resnum = int(parts[4])
                                if r465_chain == chain_id:
                                    missing_resnums.add(r465_resnum)
                            except (ValueError, IndexError):
                                pass

            # --- 3. Build full sequence and mask from SEQRES
            seqres_residues = seqres_by_chain.get(chain_id, [])
            if not seqres_residues:
                return "", [], chain_id

            # Build the SEQRES-to-PDB-resnum mapping.
            # When REMARK 465 data is available, combine resolved + missing
            # residue numbers to get the correct mapping (handles negative
            # numbering and non-standard offsets like 9JFS chain B: -52..193).
            all_resnums = sorted(resolved_resnums | missing_resnums)

            full_sequence_chars: List[str] = []
            resolved_mask: List[bool] = []

            for i, resname in enumerate(seqres_residues):
                one_letter = _resname_to_one_letter(resname)
                full_sequence_chars.append(one_letter)
                if all_resnums and len(all_resnums) == len(seqres_residues):
                    # Use proper PDB residue number mapping
                    resolved_mask.append(all_resnums[i] in resolved_resnums)
                else:
                    # Fallback: assume 1-based SEQRES == PDB resnum
                    resolved_mask.append((i + 1) in resolved_resnums)

            return "".join(full_sequence_chars), resolved_mask, chain_id

    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error extracting full sequence from PDB: {e}")
        return "", [], ""


def get_full_sequence_from_cif(
    cif_file: str, chain_id: Optional[str] = None
) -> Tuple[str, List[bool], str]:
    """
    Extract the full deposited sequence and resolved-residue mask from mmCIF.

    Reads ``_entity_poly.pdbx_seq_one_letter_code_can`` for the complete
    sequence and ``_pdbx_unobs_or_zero_occ_residues`` for missing positions.

    Args:
        cif_file: Path to mmCIF file (may be gzip-compressed).
        chain_id: Auth chain ID to extract.  If ``None`` the first RNA chain
            is used.

    Returns:
        Same as :func:`get_full_sequence_from_pdb`.
    """
    try:
        import re as _re  # pylint: disable=import-outside-toplevel

        with open_structure_file(Path(cif_file), "r") as fh:
            cif_text = fh.read()

        # --- 1. Find entity_poly entries and pick the RNA entity
        #     Look for _entity_poly.pdbx_seq_one_letter_code_can
        #     and _entity_poly.pdbx_strand_id to map entity → chains
        entity_blocks = _re.findall(
            r"_entity_poly\.pdbx_strand_id\s+(\S+)",
            cif_text,
        )
        entity_seqs = _re.findall(
            r"_entity_poly\.pdbx_seq_one_letter_code_can\s*\n?;([^;]+);",
            cif_text,
        )
        entity_types = _re.findall(
            r"_entity_poly\.type\s+['\"]?([^'\";\n]+)",
            cif_text,
        )

        if not entity_seqs:
            # Try single-line format
            entity_seqs = _re.findall(
                r"_entity_poly\.pdbx_seq_one_letter_code_can\s+(\S+)",
                cif_text,
            )

        # Match entities to types and chains
        target_seq = ""
        target_chains = ""

        for i, etype in enumerate(entity_types):
            if "ribonucleotide" in etype.lower():
                if i < len(entity_seqs):
                    target_seq = entity_seqs[i].replace("\n", "").replace(" ", "")
                if i < len(entity_blocks):
                    target_chains = entity_blocks[i]
                break

        if not target_seq:
            return "", [], ""

        # Resolve chain
        available_chains = [c.strip() for c in target_chains.split(",")]
        if chain_id is None:
            chain_id = available_chains[0] if available_chains else ""
        if not chain_id:
            return "", [], ""

        # --- 2. Find unobserved residues from _pdbx_unobs_or_zero_occ_residues
        #     This is a loop category; parse the relevant columns.
        unobs_positions: Set[int] = set()

        # Try to find the loop block
        unobs_match = _re.search(
            r"loop_\s*\n"
            r"((?:_pdbx_unobs_or_zero_occ_residues\.\w+\s*\n)+)"
            r"((?:(?!loop_|#|_).*\n)*)",
            cif_text,
        )
        if unobs_match:
            header_text = unobs_match.group(1)
            data_text = unobs_match.group(2)

            headers = _re.findall(
                r"_pdbx_unobs_or_zero_occ_residues\.(\w+)", header_text
            )
            auth_asym_idx = (
                headers.index("auth_asym_id") if "auth_asym_id" in headers else None
            )
            seq_id_idx = (
                headers.index("auth_seq_id") if "auth_seq_id" in headers else None
            )
            label_seq_idx = (
                headers.index("label_seq_id") if "label_seq_id" in headers else None
            )

            # Use label_seq_id (entity sequence numbering) if available
            pos_idx = label_seq_idx if label_seq_idx is not None else seq_id_idx

            for data_line in data_text.strip().split("\n"):
                fields = data_line.split()
                if not fields or len(fields) <= max(
                    f for f in [auth_asym_idx, pos_idx] if f is not None
                ):
                    continue
                if auth_asym_idx is not None and fields[auth_asym_idx] != chain_id:
                    continue
                if pos_idx is not None:
                    try:
                        unobs_positions.add(int(fields[pos_idx]))
                    except ValueError:
                        pass

        # --- 3. Build mask
        full_sequence = target_seq
        resolved_mask = [
            (i + 1) not in unobs_positions for i in range(len(full_sequence))
        ]

        return full_sequence, resolved_mask, chain_id

    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error extracting full sequence from CIF: {e}")
        return "", [], ""


def get_full_sequence(
    structure_file: str, chain_id: Optional[str] = None
) -> Tuple[str, List[bool], str]:
    """
    Extract full deposited sequence and resolved mask from any structure file.

    Dispatches to PDB or CIF parser based on file extension.

    Args:
        structure_file: Path to PDB or mmCIF file.
        chain_id: Optional chain ID.

    Returns:
        Tuple of ``(full_sequence, resolved_mask, used_chain_id)``.
    """
    path = Path(structure_file)
    # Strip .gz to check real extension
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]

    if name.endswith(".cif"):
        return get_full_sequence_from_cif(structure_file, chain_id)
    return get_full_sequence_from_pdb(structure_file, chain_id)


def remap_dot_bracket(
    resolved_dot_bracket: str,
    resolved_mask: List[bool],
) -> str:
    """
    Expand a dot-bracket string from resolved-only indexing to full-sequence indexing.

    Inserts ``'.'`` at every unresolved position so the result has the same
    length as the full deposited sequence.

    Args:
        resolved_dot_bracket: Dot-bracket from base-pair extraction (covers
            only resolved nucleotides).
        resolved_mask: Boolean list for the full sequence (``True`` = resolved).

    Returns:
        Expanded dot-bracket string with ``len(resolved_mask)`` characters.

    Raises:
        ValueError: If the number of ``True`` values in *resolved_mask*
            doesn't match *len(resolved_dot_bracket)*.
    """
    resolved_count = sum(resolved_mask)
    if resolved_count != len(resolved_dot_bracket):
        raise ValueError(
            f"Resolved count ({resolved_count}) != "
            f"dot-bracket length ({len(resolved_dot_bracket)})"
        )

    result = []
    res_idx = 0
    for is_resolved in resolved_mask:
        if is_resolved:
            result.append(resolved_dot_bracket[res_idx])
            res_idx += 1
        else:
            result.append(".")
    return "".join(result)


# pylint: disable=too-many-branches, too-many-statements
def parse_fr3d_basepairs(
    basepair_file: str,
    unit_id_to_position: Dict[str, int],
    sequence_length: int,
    include_pseudoknots: bool = False,
) -> str:
    """
    Parse FR3D basepair output and convert to dot-bracket notation.

    FR3D outputs lines like:
    1S72|1|0|A|1193    cWW    1S72|1|0|U|1205    0

    We consider Watson-Crick and wobble pairs (cWW). Nested pairs (crossing=0)
    get () brackets. Pseudoknots (crossing>0) get Aa, Bb, Cc notation.

    Args:
        basepair_file: Path to FR3D basepair output file.
        unit_id_to_position: Mapping from unit_id to sequence position (0-based).
        sequence_length: Length of the sequence.
        include_pseudoknots: If True, include pseudoknots using Aa, Bb notation.

    Returns:
        Dot-bracket notation string (with extended alphabet if pseudoknots included).
    """
    if sequence_length == 0:
        return ""

    dot_bracket = ["."] * sequence_length

    # Collect all base pairs: (pos1, pos2, crossing_number)
    all_pairs: List[Tuple[int, int, int]] = []
    seen_pairs: Set[Tuple[int, int]] = set()

    with open(basepair_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 4:
                continue

            unit_id1, interaction, unit_id2, crossing_str = (
                parts[0],
                parts[1],
                parts[2],
                parts[3],
            )

            # Only consider Watson-Crick pairs (cWW)
            if interaction not in ["cWW", "cWw", "cwW"]:
                continue

            try:
                crossing = int(crossing_str)
            except ValueError:
                continue

            # Skip pseudoknots if not requested
            if crossing != 0 and not include_pseudoknots:
                continue

            # Get sequence positions
            pos1 = _get_position_from_unit_id(unit_id1, unit_id_to_position)
            pos2 = _get_position_from_unit_id(unit_id2, unit_id_to_position)

            if pos1 is None or pos2 is None:
                continue

            if pos1 == pos2:
                continue

            # Normalize pair order and avoid duplicates
            pair = (min(pos1, pos2), max(pos1, pos2))
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)

            all_pairs.append((pair[0], pair[1], crossing))

    # Separate nested pairs (crossing=0) from pseudoknots (crossing>0)
    nested_pairs = [(p1, p2) for p1, p2, c in all_pairs if c == 0]
    pseudoknot_pairs = [(p1, p2, c) for p1, p2, c in all_pairs if c > 0]

    # Track which positions are already paired
    paired_positions: Set[int] = set()

    # First assign () to nested pairs
    for pos1, pos2 in nested_pairs:
        if pos1 in paired_positions or pos2 in paired_positions:
            continue
        dot_bracket[pos1] = "("
        dot_bracket[pos2] = ")"
        paired_positions.add(pos1)
        paired_positions.add(pos2)

    # Then assign pseudoknots using Aa, Bb, Cc notation
    if include_pseudoknots and pseudoknot_pairs:
        # Group pseudoknots by crossing number to assign different letters
        pk_by_crossing: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
        for pos1, pos2, crossing in pseudoknot_pairs:
            pk_by_crossing[crossing].append((pos1, pos2))

        # Assign letters: A/a for first group, B/b for second, etc.
        pk_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        for idx, (crossing, pairs) in enumerate(sorted(pk_by_crossing.items())):
            if idx >= len(pk_letters):
                break  # Only support 26 pseudoknot levels

            open_char = pk_letters[idx]
            close_char = pk_letters[idx].lower()

            for pos1, pos2 in pairs:
                if pos1 in paired_positions or pos2 in paired_positions:
                    continue
                dot_bracket[pos1] = open_char
                dot_bracket[pos2] = close_char
                paired_positions.add(pos1)
                paired_positions.add(pos2)

    result = "".join(dot_bracket)

    # Clean up invalid patterns like "()" - only for regular brackets
    result = _fix_invalid_pairs(result)

    return result


def _get_position_from_unit_id(
    unit_id: str,
    unit_id_to_position: Dict[str, int],
) -> Optional[int]:
    """
    Get sequence position from a unit_id.

    First tries exact match, then tries matching by key components.

    Args:
        unit_id: FR3D unit identifier.
        unit_id_to_position: Mapping dictionary.

    Returns:
        Sequence position (0-based) or None if not found.
    """
    # Try exact match first
    if unit_id in unit_id_to_position:
        return unit_id_to_position[unit_id]

    # Try matching by PDB|model|chain|resname|resnum pattern
    parts = unit_id.split("|")
    if len(parts) >= 5:
        # Try without model/symmetry variations
        for key in unit_id_to_position:
            key_parts = key.split("|")
            if len(key_parts) >= 5:
                # Match by chain and residue number
                if parts[2] == key_parts[2] and parts[4] == key_parts[4]:
                    return unit_id_to_position[key]

    return None


def _fix_invalid_pairs(dot_bracket: str) -> str:
    """
    Fix invalid dot-bracket patterns like "()" by converting them to "..".

    Args:
        dot_bracket: Dot-bracket notation string.

    Returns:
        Cleaned dot-bracket notation.
    """
    result = list(dot_bracket)
    changed = True

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


# pylint: disable=too-many-arguments,too-many-positional-arguments
def get_secondary_structure_fr3d(
    structure_file: str,
    output_dir: str,
    chain_id: Optional[str] = None,
    model_id: Optional[int] = None,
    include_pseudoknots: bool = False,
    quiet: bool = False,
) -> Tuple[Optional[str], Optional[str]]:
    """
    Extract RNA secondary structure from a 3D structure using FR3D.

    This is the main entry point for FR3D-based structure extraction.

    Args:
        structure_file: Path to the structure file (.cif or .pdb).
        output_dir: Directory for intermediate files.
        chain_id: Optional chain ID to extract.
        model_id: Optional model ID to extract. If None, uses first model.
        include_pseudoknots: If True, include pseudoknots using Aa, Bb notation.
        quiet: If True, suppress verbose output.

    Returns:
        Tuple of (sequence, dot_bracket) or (None, None) if failed.
    """
    structure_path = Path(structure_file)

    # Extract sequence based on file type
    if structure_path.suffix.lower() == ".cif":
        sequence, unit_id_to_position = extract_sequence_from_cif(
            str(structure_path), chain_id, model_id, quiet=quiet
        )
    else:
        sequence, unit_id_to_position = extract_sequence_from_pdb(
            str(structure_path), chain_id, quiet=quiet
        )

    if not sequence:
        if not quiet:
            print(f"Could not extract sequence from {structure_file}")
        return None, None

    if not quiet:
        print(f"Extracted sequence: {len(sequence)} nucleotides")

    # Run FR3D to get basepair annotations
    basepair_file = run_fr3d(str(structure_path), output_dir, quiet=quiet)
    if not basepair_file:
        if not quiet:
            print("FR3D analysis failed, returning unpaired structure")
        return sequence, "." * len(sequence)

    # Parse basepairs and convert to dot-bracket
    dot_bracket = parse_fr3d_basepairs(
        basepair_file, unit_id_to_position, len(sequence), include_pseudoknots
    )

    pk_count = sum(1 for c in dot_bracket if c.isupper() and c not in ".()")
    if not quiet:
        print(
            f"Secondary structure: {dot_bracket.count('(')} base pairs, {pk_count} pseudoknot pairs"
        )
    return sequence, dot_bracket
