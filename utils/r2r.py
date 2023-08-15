"""
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

import os

R2R_META = "r2r.meta.txt"
R2R_INPUT = "r2r-input.sto.txt"
R2R_GSC = "r2r-input.gsc.sto.txt"
R2R_SVG = "r2r.svg"
R2R_SEQ_ID = "input"


def run_traveler(fasta_input, output_folder, seq_id):
    """Run traveler."""
    traveler_params = (
        f"--template-structure --file-format traveler "
        f"{output_folder}/traveler-template.xml "
        f"{fasta_input}"
    )
    cmd = (
        "traveler --verbose "
        f"--target-structure {fasta_input} {traveler_params} "
        f"--all {output_folder}/{seq_id} > {output_folder}/traveler.log"
    )
    os.system(cmd)


def parse_fasta(fasta_input):
    """Parse a fasta file and return the sequence and structure."""
    with open(fasta_input, "r", encoding="utf-8") as f_in:
        lines = f_in.readlines()
    sequence = ""
    structure = ""
    for i, line in enumerate(lines):
        if i == 0 and line.startswith(">"):
            seq_id = line.replace(">", "").split(" ")[0].strip()
        if i == 1:
            sequence = line.strip()
        if i == 2:
            structure = line.strip()
    if not sequence and not structure:
        raise ValueError("Invalid fasta file.")
    return seq_id, sequence, structure


def generate_r2r_input_file(sequence, structure, output_folder):
    """Generate an input file for R2R."""
    with open(os.path.join(output_folder, R2R_INPUT), "w", encoding="utf-8") as f_out:
        f_out.write("# STOCKHOLM 1.0\n")
        f_out.write(f"{R2R_SEQ_ID}          {sequence}\n")
        f_out.write(f"#=GC SS_cons   {structure}\n")
        f_out.write("#=GF R2R tick_label_disable_default_numbering\n")
        f_out.write("//\n")
    with open(os.path.join(output_folder, R2R_META), "w", encoding="utf-8") as f_out:
        f_out.write(f"{output_folder}/{R2R_GSC}\toneseq\t{R2R_SEQ_ID}\n")


def run_r2r(output_folder):
    """Run R2R: generate GSC alignment and SVG output."""
    cmd = (
        f"r2r --GSC-weighted-consensus {output_folder}/{R2R_INPUT} "
        f"{output_folder}/{R2R_GSC} 3 0.97 0.9 0.75 4 0.97 0.9 0.75 0.5 0.1"
    )
    os.system(cmd)
    cmd = (
        f"r2r --disable-usage-warning {output_folder}/{R2R_META} "
        f"{output_folder}/{R2R_SVG} > {output_folder}/r2r.log"
    )
    os.system(cmd)
    return f"{output_folder}/{R2R_SVG}"


def clean_r2r_output(output_folder):
    """Remove filename from R2R output."""
    r2r_svg_file = os.path.join(output_folder, R2R_SVG)
    with open(r2r_svg_file, "r", encoding="utf-8") as f_in:
        lines = f_in.readlines()
    with open(r2r_svg_file, "w", encoding="utf-8") as f_out:
        for line in lines:
            if R2R_GSC in line:
                continue
            f_out.write(line)
