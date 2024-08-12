"""
This module contains functions to visualise consensus
secondary structures for Rfam families using RNArtist and R-scape
layouts and compare the number of overlaps between the two.
"""

import shutil
import tempfile
from pathlib import Path

from .rfam import RfamSeed
from .runner import runner


def get_overlaps(output_folder, rfam_acc):
    """Count the number of overlaps in a Traveler XML file."""
    svg_file = Path(output_folder) / rfam_acc / f"query-{rfam_acc}.colored.svg"
    overlaps_file = Path(output_folder) / rfam_acc / f"query-{rfam_acc}.overlaps"
    if not svg_file.exists() or not overlaps_file.exists():
        print(f"Error: overlaps or svg file not found in {output_folder}")
        return float("inf")
    with open(overlaps_file, "r") as f_in:
        overlaps = int(f_in.read().strip())
    return overlaps


def get_rfam_consensus(rfam_acc):
    """Get the consensus sequence of an Rfam family."""
    rfam_seed = RfamSeed().get_rfam_seed(rfam_acc)
    consensus_seq = ""
    rfam_author_encoding = "ISO-8859-1"
    with open(rfam_seed, "r", encoding=rfam_author_encoding) as f_in:
        for line in f_in:
            if line.startswith("#=GC RF"):
                consensus_seq += line.split()[2]
    consensus_seq = consensus_seq.replace("-", "").replace(".", "").upper()
    return consensus_seq


def compare_rnartist_and_rscape(rfam_acc):
    """Run Traveler on both R-scape and RNArtist consensus 2D layouts and count overlaps."""
    # pylint: disable=consider-using-with
    query_fasta_file = Path(
        tempfile.NamedTemporaryFile(suffix=".fasta", delete=False).name
    )
    with open(query_fasta_file, "w") as f_out:
        f_out.write(">query\n")
        f_out.write(f"{get_rfam_consensus(rfam_acc)}\n")

    rscape_output = tempfile.mkdtemp()
    rscape_cmd = f"r2dt.py rfam draw {rfam_acc} {query_fasta_file} {rscape_output} --quiet --rscape"
    runner.run(rscape_cmd)
    rscape_overlaps = get_overlaps(rscape_output, rfam_acc)
    shutil.rmtree(rscape_output)

    if rscape_overlaps > 0:
        rnartist_output = tempfile.mkdtemp()
        rnartist_cmd = (
            f"r2dt.py rfam draw {rfam_acc} {query_fasta_file} "
            f"{rnartist_output} --rnartist --quiet"
        )
        runner.run(rnartist_cmd, print_output=True)
        rnartist_overlaps = get_overlaps(rnartist_output, rfam_acc)
        shutil.rmtree(rnartist_output)
    else:
        rnartist_overlaps = rscape_overlaps

    query_fasta_file.unlink()
    query_fasta_file.with_suffix(".fasta.ssi").unlink(missing_ok=True)

    if rnartist_overlaps < rscape_overlaps:
        return "rnartist", {"rnartist": rnartist_overlaps, "rscape": rscape_overlaps}
    return "rscape", {}
