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
import re

import requests  # pylint: disable=import-error
import RNA  # pylint: disable=import-error
from colorhash import ColorHash  # pylint: disable=import-error

MAX_INSERTIONS = 100


def get_r2dt_version_header():
    """Return a welcome message including release information."""
    header = """# R2DT :: visualise RNA secondary structure using templates
# Version 1.3 (October 2022)
# https://github.com/RNAcentral/R2DT
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"""
    return header


def get_ribotyper_output(fasta_input, output_folder, cm_library, skip_ribovore_filters):
    """
    Run ribotyper on the fasta sequences to select the best matching covariance
    model.
    """
    ribotyper_long_out = os.path.join(
        output_folder, os.path.basename(output_folder) + ".ribotyper.long.out"
    )
    if not os.path.exists(ribotyper_long_out):
        cmd = (
            f"ribotyper.pl --skipval -i {cm_library}/modelinfo.txt "
            f"-f {fasta_input} {output_folder}"
        )
        print(cmd)
        os.system(cmd)
    f_out = os.path.join(output_folder, "hits.txt")
    if not skip_ribovore_filters:
        cmd = (
            f"cat {ribotyper_long_out} | grep -v '^#' | "
            f"grep -v MultipleHits | grep PASS | "
            f"awk -v OFS='\t' '{{print $2, $8, $3}}' > {f_out}"
        )
    else:
        cmd = (
            f"cat {ribotyper_long_out} | grep -v '^#' "
            f"| grep -v NoHits | "
            f"awk -v OFS='\t' '{{print $2, $8, $3}}' > {f_out}"
        )
    os.system(cmd)
    return f_out


def remove_large_insertions_pfam_stk(filename):
    """
    The Pfam Stockholm files can contain 9 or 11 lines depending on whether
    the description line is present.
    """
    with open(filename, "r", encoding="utf-8") as f_stockholm:
        lines = f_stockholm.readlines()
        if len(lines) == 9:
            sequence = lines[3]
            gr_pp = lines[4]
            gr_ss = lines[5]
            gc_ss_cons = lines[6]
            gc_rf = lines[7]
        elif len(lines) == 11:
            sequence = lines[5]
            gr_pp = lines[6]
            gr_ss = lines[7]
            gc_ss_cons = lines[8]
            gc_rf = lines[9]
        elif len(lines) == 3:
            sequence = lines[0]
            gr_pp = ""
            gr_ss = lines[1]
            gc_ss_cons = lines[2]
            gc_rf = ""
        else:
            print("Unexpected number of lines in pfam stk")
            return
        # the tilda and period characters represent insert states in WUSS notation
        match = re.finditer(r"([\.~_:]{" + str(MAX_INSERTIONS) + ",})", gc_ss_cons)
        if not match:
            return
        for span in match:
            sequence = (
                sequence[: span.start()]
                + "@" * (span.end() - span.start())
                + sequence[span.end() :]
            )
            if gr_pp:
                gr_pp = (
                    gr_pp[: span.start()]
                    + "@" * (span.end() - span.start())
                    + gr_pp[span.end() :]
                )
            gr_ss = (
                gr_ss[: span.start()]
                + "@" * (span.end() - span.start())
                + gr_ss[span.end() :]
            )
            gc_ss_cons = (
                gc_ss_cons[: span.start()]
                + "@" * (span.end() - span.start())
                + gc_ss_cons[span.end() :]
            )
            if gc_rf:
                gc_rf = (
                    gc_rf[: span.start()]
                    + "@" * (span.end() - span.start())
                    + gc_rf[span.end() :]
                )
        if len(lines) == 9:
            lines[3] = re.sub(r"@+", "XXXX", sequence)
            lines[4] = re.sub(r"@+", "XXXX", gr_pp)
            lines[5] = re.sub(r"@+", "~~~~", gr_ss)
            lines[6] = re.sub(r"@+", "~~~~", gc_ss_cons)
            lines[7] = re.sub(r"@+", "xxxx", gc_rf)
        elif len(lines) == 11:
            lines[5] = re.sub(r"@+", "XXXX", sequence)
            lines[6] = re.sub(r"@+", "XXXX", gr_pp)
            lines[7] = re.sub(r"@+", "~~~~", gr_ss)
            lines[8] = re.sub(r"@+", "~~~~", gc_ss_cons)
            lines[9] = re.sub(r"@+", "xxxx", gc_rf)
        elif len(lines) == 3:
            lines[0] = re.sub(r"@+", "XXXX", sequence)
            lines[1] = re.sub(r"@+", "~~~~", gr_ss)
            lines[2] = re.sub(r"@+", "~~~~", gc_ss_cons)
    with open(filename, "w", encoding="utf-8") as f_stockholm:
        for line in lines:
            f_stockholm.write(line)


def get_insertions(filename):
    """Extract insertions to be folded."""
    sequence = ""
    with open(filename, "r", encoding="utf-8") as f_stockholm:
        lines = f_stockholm.readlines()
        if len(lines) == 9:
            sequence = lines[3].split()[1]
        elif len(lines) == 11:
            sequence = lines[5].split()[1]
        trimmed = sequence.replace("-", "")
    match = re.finditer(r"[a-z]+", trimmed)
    return match


# pylint: disable-next=too-many-locals
def get_full_constraint(filename):
    """Get folding constraint from an Infernal alignment."""
    constraint = ""
    with open(filename, "r", encoding="utf-8") as f_stockholm:
        lines = f_stockholm.readlines()
        if len(lines) == 9:
            gc_ss = lines[5].split()[3]
            sequence = lines[3].split()[1]
        elif len(lines) == 11:
            gc_ss = lines[7].split()[3]
            sequence = lines[5].split()[1]
    deletions = re.finditer("-+", sequence)
    sub_ss = list(gc_ss)
    trimmed = sequence.replace("-", "")
    match = re.finditer(r"[a-z]+", trimmed)
    adjust = 0
    for span in deletions:
        i = span.start() - adjust
        j = span.end() - adjust
        sub_ss = sub_ss[0:i] + sub_ss[j : len(sub_ss)]
        adjust = adjust + j - i
    match_arr = set()
    for span in match:
        i = span.start()
        j = span.end()
        while i < j:
            match_arr.add(i)
            i += 1
    for i, val in enumerate(sub_ss):
        if i in match_arr:
            constraint += "."
        elif val in ["<", "(", "{"]:
            constraint += "("
        elif val in [">", ")", "}"]:
            constraint += ")"
        elif val in [",", "-", ":"]:
            constraint += "x"
        else:
            constraint += "x"
    return constraint


# pylint: disable-next=too-many-locals
def fold_insertions_only(sequence, constraint, filename):
    """Use RNAfold to fold insertions."""
    match = get_insertions(filename)
    list_seq = list(sequence)
    list_con = list(constraint)
    final_list = ""
    for span in match:
        i = span.start()
        j = span.end()
        while i > 1 and (list_con[i - 1] == "x" or list_con[i - 1] == "."):
            i -= 1
        while j < (len(list_con) - 2) and (list_con[j] == "x" or list_con[j] == "."):
            j += 1
        subsequence = "".join(list_seq[i:j])
        constart = list_con[0:i]
        subconstraint = "".join(list_con[i:j])
        formatted_constraint = subconstraint.replace("(", "x").replace(")", "x")
        conend = list_con[j : len(list_con) + 1]
        model_details = RNA.md()
        model_details.min_loop_size = 0
        fold_compound = RNA.fold_compound(subsequence, model_details)
        constraint_options = (
            RNA.CONSTRAINT_DB | RNA.CONSTRAINT_DB_ENFORCE_BP | RNA.CONSTRAINT_DB_DEFAULT
        )
        fold_compound.hc_add_from_db(formatted_constraint, constraint_options)
        (secondary_structure, min_free_energy) = fold_compound.mfe()
        if min_free_energy < 99999:
            list_con = constart
            for i, _ in enumerate(secondary_structure):
                if subconstraint[i] != ".":
                    list_con += subconstraint[i]
                else:
                    list_con += secondary_structure[i]
            list_con += conend
        else:
            print("Structure exceeds limits of RNAfold, ignoring constraint")
            return None
    final_list = "".join(list_con).replace("x", ".")
    return final_list


def handle_exclusion(exclusion, r2dt_constraint):
    """Process secondary structure exclusion file."""
    try:
        with open(exclusion, "r", encoding="utf-8") as f_exclusion:
            exclusion_string = f_exclusion.read().strip()
            if len(exclusion_string) != len(r2dt_constraint):
                print("Exclusion ignored, not same length as sequence")
                return r2dt_constraint
            if re.search(r"[^\.x]", exclusion_string):
                print("Invalid characters in exclusion string, should only be . and x")
                return r2dt_constraint
            temp_constraint = ""
            for i, _ in enumerate(exclusion_string):
                if exclusion_string[i] == "x":
                    if r2dt_constraint[i] != "." and r2dt_constraint[i] != "x":
                        print(
                            "Invalid exclusion in position "
                            + str(i)
                            + " conflicts with template, constraint ignored"
                        )
                        temp_constraint += r2dt_constraint[i]
                    else:
                        temp_constraint += "x"
                else:
                    temp_constraint += r2dt_constraint[i]
            return temp_constraint
    except FileNotFoundError:
        print("Constraint file not found")
        return r2dt_constraint


# pylint: disable=too-many-arguments, too-many-locals, too-many-branches, too-many-statements
def fold_insertions(input_fasta, exclusion, source, filename, model_id, fold_type):
    """Fold insertions using RNAfold."""
    with open(input_fasta, "r", encoding="utf-8") as f_fasta:
        orig = f_fasta.readlines()
    r2dt_constraint = orig[2].strip()
    orig_sequence = orig[1].strip()
    secondary_structure = ""
    if fold_type not in [
        "insertions_only",
        "full_molecule",
        "all_constraints_enforced",
    ]:
        if source == "rfam":
            insertions_only = [
                "rRNA",
                "lncRNA",
                "Cis-regulatory element",
                "ribozyme",
                "CRISPR",
                "antisense",
                "antitoxin",
                "Intron",
            ]
            full_molecule = ["snRNA", "snoRNA", "sRNA", "tRNA", "miRNA"]
            rfam_api = requests.get(
                f"http://rfam.org/family/{model_id}?content-type=application/json",
                timeout=60,
            )
            if any(
                x in rfam_api.json()["rfam"]["curation"]["type"] for x in full_molecule
            ):
                fold_type = "full_molecule"
            elif any(
                x in rfam_api.json()["rfam"]["curation"]["type"]
                for x in insertions_only
            ):
                fold_type = "insertions_only"
            else:
                fold_type = "insertions_only"
        elif source == "gtrnadb":
            fold_type = "full_molecule"
        else:
            fold_type = "insertions_only"
    if fold_type == "all_constraints_enforced":
        constraint = get_full_constraint(filename)
        final_constraint = ""
        if exclusion:
            final_constraint = handle_exclusion(exclusion, constraint)
        else:
            final_constraint = constraint
        model_details = RNA.md()
        model_details.min_loop_size = 0
        fold_compound = RNA.fold_compound(orig_sequence, model_details)
        constraint_options = (
            RNA.CONSTRAINT_DB | RNA.CONSTRAINT_DB_ENFORCE_BP | RNA.CONSTRAINT_DB_DEFAULT
        )
        fold_compound.hc_add_from_db(final_constraint, constraint_options)
        (secondary_structure, min_free_energy) = fold_compound.mfe()
        if min_free_energy > 99999:
            print("Structure exceeds limits of RNAfold, ignoring constraint")
            return
    elif fold_type == "full_molecule":
        if exclusion:
            final_constraint = handle_exclusion(exclusion, r2dt_constraint)
        else:
            final_constraint = r2dt_constraint
        model_details = RNA.md()
        model_details.min_loop_size = 0
        fold_compound = RNA.fold_compound(orig_sequence, model_details)
        constraint_options = (
            RNA.CONSTRAINT_DB | RNA.CONSTRAINT_DB_ENFORCE_BP | RNA.CONSTRAINT_DB_DEFAULT
        )
        fold_compound.hc_add_from_db(final_constraint, constraint_options)
        (secondary_structure, min_free_energy) = fold_compound.mfe()
        if min_free_energy > 99999:
            print("Structure exceeds limits of RNAfold, ignoring constraint")
            return
    # Default option
    elif fold_type == "insertions_only":
        constraint = ""
        if exclusion:
            constraint = handle_exclusion(exclusion, r2dt_constraint)
        else:
            constraint = r2dt_constraint
        secondary_structure = fold_insertions_only(orig_sequence, constraint, filename)
    if secondary_structure:
        with open(input_fasta, "w", encoding="utf-8") as f_fasta:
            f_fasta.write(orig[0])
        constraint_differences = ""
        for i, val in enumerate(r2dt_constraint):
            if val == secondary_structure[i]:
                constraint_differences += "-"
            else:
                constraint_differences += "*"
        with open(input_fasta, "a", encoding="utf-8") as f_fasta:
            f_fasta.write(orig_sequence)
            f_fasta.write("\n" + secondary_structure)
            f_fasta.write("\n" + constraint_differences)


def get_infernal_posterior_probabilities(input_file, output_file):
    """
    Parse an Infernal stockholm alignment in Pfam format.
    Return a tsv file with posterior probabilities
    that can be used to propagate the values to JSON
    and colour the SVGs.
    """
    sequence = ""
    post_prob = ""
    with open(input_file, "r", encoding="utf-8") as f_in:
        for line in f_in:
            match = re.match(r"^#=GR\s+.+?\s+PP\s+([123456789*.]+)$", line)
            if match:
                post_prob = match.group(1)
                continue
            if line.startswith("#="):
                continue
            match = re.match(r"^.+?\s{2,}(.+?)$", line)
            if match:
                sequence = match.group(1)
                continue
    with open(output_file, "w", encoding="utf-8") as f_out:
        header = ["residue_index", "residue_name", "posterior_probability"]
        f_out.write("\t".join(header) + "\n")
        index = 1
        for nucleotide, prob in zip(list(sequence), list(post_prob)):
            if nucleotide == "-":
                continue
            if prob == "*":
                prob = 10
            f_out.write(f"{index}\t{nucleotide}\t{prob}\n")
            index += 1


def generate_thumbnail(image, description):
    """Generate a thumbnail SVG as an outline of the 2D diagram."""
    move_to_start_position = None
    color = ColorHash(description).hex
    points = []
    for _, line in enumerate(image.split("\n")):
        if "width" in line and not "stroke-width" in line:
            width = re.findall(r'width="(\d+(\.\d+)?)"', line)
        if "height" in line:
            height = re.findall(r'height="(\d+(\.\d+)?)"', line)
        for nt_block in re.finditer(
            r'<text x="(\d+)(\.\d+)?" y="(\d+)(\.\d+)?".*?</text>', line
        ):
            if "numbering-label" in nt_block.group(0):
                continue
            if not move_to_start_position:
                move_to_start_position = f"M{nt_block.group(1)} {nt_block.group(3)} "
            points.append(f"L{nt_block.group(1)} {nt_block.group(3)}")
    if len(points) < 200:
        stroke_width = "3"
    elif len(points) < 500:
        stroke_width = "4"
    elif len(points) < 3000:
        stroke_width = "4"
    else:
        stroke_width = "2"
    thumbnail = (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{width[0][0]}" height="{height[0][0]}">'
        f'<path style="stroke:{color};stroke-width:{stroke_width}px;'
        f'fill:none;" d="'
    )
    thumbnail += move_to_start_position
    thumbnail += " ".join(points)
    thumbnail += '"/></svg>'
    return thumbnail
