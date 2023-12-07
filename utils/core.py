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
from pathlib import Path

from rich import print as rprint

from . import config, gtrnadb, rfam, shared
from .runner import runner


# pylint: disable=too-many-arguments
# pylint: disable=too-many-branches
# pylint: disable=too-many-locals
# pylint: disable=too-many-statements
# pylint: disable=too-many-return-statements
def visualise(
    rna_type,
    fasta_input,
    output_folder,
    seq_id,
    model_id,
    constraint,
    exclusion,
    fold_type,
    domain=None,
    isotype=None,
    start=None,
    end=None,
):
    """Main visualisation routine that invokes Traveler."""
    if model_id:
        rprint(f"Visualising {seq_id} with {model_id}")
    else:
        rprint(f"Visualising {seq_id} with {domain} {isotype}")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    filename_template = os.path.join(output_folder, f"{seq_id}_type.txt")
    if rna_type.lower() == "lsu":
        cm_library = config.RIBOVISION_LSU_CM_LIBRARY
        template_layout = config.RIBOVISION_LSU_TRAVELER
        template_structure = config.RIBOVISION_LSU_BPSEQ
    elif rna_type.lower() == "ssu":
        cm_library = config.RIBOVISION_SSU_CM_LIBRARY
        template_layout = config.RIBOVISION_SSU_TRAVELER
        template_structure = config.RIBOVISION_SSU_BPSEQ
    elif rna_type.lower() == "rnasep":
        cm_library = config.RNASEP_CM_LIBRARY
        template_layout = config.RNASEP_TRAVELER
        template_structure = config.RNASEP_BPSEQ
    elif rna_type.lower() == "crw":
        cm_library = config.CRW_CM_LIBRARY
        template_layout = config.CRW_PS_LIBRARY
        template_structure = config.CRW_FASTA_LIBRARY
    elif rna_type.lower() == "rfam":
        if not model_id.startswith("RF"):
            model_id = rfam.get_rfam_acc_by_id(model_id)
        temp_sto_unfiltered = filename_template.replace("type", "unfiltered")
        temp_acc_list = filename_template.replace("type", "seed_list")
    elif rna_type.lower() == "local_data":
        cm_library = os.path.join(config.LOCAL_DATA, model_id)
        template_layout = cm_library
        template_structure = cm_library
    elif rna_type.lower() == "gtrnadb":
        model_id = domain + "_" + isotype
    else:
        rprint("Please specify RNA type")
        return

    temp_fasta = filename_template.replace("type", "fasta")
    temp_sto = filename_template.replace("type", "sto")
    temp_depaired = filename_template.replace("type", "depaired")
    temp_stk = filename_template.replace("type", "stk")
    temp_stk_original = filename_template.replace("type", "stk_original")
    temp_post_prob = filename_template.replace("type", "post_prob")
    temp_pfam_stk = filename_template.replace("type", "pfam_stk")
    temp_pfam_stk_original = filename_template.replace("type", "pfam_stk_original")
    temp_afa = filename_template.replace("type", "afa")
    temp_afa_original = filename_template.replace("type", "afa_original")
    temp_map = filename_template.replace("type", "map")

    # get sequence from fasta file
    seq_range = f"-c {start}..{end}" if start and end else ""
    if not os.path.exists(f"{fasta_input}.ssi"):
        runner.run(f"esl-sfetch --index {fasta_input}")
    cmd = f"esl-sfetch {seq_range} {fasta_input} {seq_id} > {temp_fasta}"
    result = runner.run(cmd)
    if result:
        raise ValueError(f"Failed esl-sfetch for: {seq_id} in {cmd}")

    # check that the model exists
    if rna_type == "rfam":
        model_path = rfam.get_rfam_cm(model_id)
        template_layout = rfam.get_traveler_template_xml(model_id)
        template_structure = rfam.get_traveler_fasta(model_id)
        # download seed alignment and list its accessions
        rfam_seed = rfam.download_rfam_seed(model_id)
        cmd = f"esl-alistat --list {temp_acc_list} {rfam_seed} > /dev/null"
        runner.run(cmd)
    elif rna_type == "gtrnadb":
        model_path = gtrnadb.get_trnascan_cm(domain, isotype)
        if not model_path:
            rprint(f"Covariance model not found for {domain} {isotype}")
            return
        template_layout = gtrnadb.get_traveler_template_xml(domain, isotype)
        template_structure = gtrnadb.get_traveler_fasta(domain, isotype)
    elif rna_type == "local_data":
        model_path = os.path.join(config.LOCAL_DATA, model_id, model_id + ".cm")
        if not os.path.exists(model_path):
            rprint(f"Model not found {model_path}")
            return
        template_layout = os.path.join(template_layout, model_id + ".xml")
        template_structure = os.path.join(template_structure, model_id + ".fasta")
    else:
        model_path = os.path.join(cm_library, model_id + ".cm")
        if not os.path.exists(model_path):
            rprint(f"Model not found {model_path}")
            return

    # align sequence to the model
    cm_options = ["", "--mxsize 2048 --maxtau 0.49"]
    for options in cm_options:
        if rna_type == "rfam":
            cmd = (
                f"cmalign --mapali {rfam_seed} --mapstr {options} "
                f"{model_path} {temp_fasta} > {temp_sto_unfiltered}"
            )
        else:
            cmd = f"cmalign {options} {model_path} {temp_fasta} > {temp_sto}"
        result = runner.run(cmd)
        if not result:
            break
    else:
        rprint(f"Failed cmalign of {seq_id} to {model_id}")
        return

    if rna_type == "rfam":
        cmd = (
            f"esl-alimanip --seq-r {temp_acc_list} {temp_sto_unfiltered} | "
            f"esl-reformat --keeprf --mingap --informat stockholm stockholm - > "
            f"{temp_sto}"
        )
        runner.run(cmd)

    # remove non-canonical Watson-Crick basepairs (e.g. C:A in URS000008DB9C_7227)
    cmd = f"esl-alidepair.pl --nc 0.5 {temp_sto} {temp_depaired}"
    result = runner.run(cmd)
    if result:
        rprint(f"Failed esl-alidepair for {seq_id}")

    has_conserved_structure = False
    with open(temp_sto) as f_stockholm:
        for line in f_stockholm.readlines():
            if line.startswith("#=GC SS_cons"):
                if "<" in line:
                    has_conserved_structure = True
                else:
                    rprint("This RNA does not have a conserved structure")
                break
    if not has_conserved_structure:
        return

    # impose consensus secondary structure and convert to pfam format
    cmd = (
        f"esl-alimanip --rna --sindi " f"--outformat pfam {temp_depaired} > {temp_stk}"
    )
    result = runner.run(cmd)
    if result:
        rprint(f"Failed esl-alimanip for {seq_id} {model_id}")
        return

    # impose consensus secondary structure and convert to pfam format
    cmd = (
        f"esl-alimanip --rna --sindi "
        f"--outformat pfam {temp_sto} > {temp_stk_original}"
    )
    result = runner.run(cmd)
    if result:
        rprint(f"Failed esl-alimanip for {seq_id} {model_id}")
        return

    # store posterior probabilities in tsv file
    shared.get_infernal_posterior_probabilities(temp_stk, temp_post_prob)

    # convert nts that are in RF-gap columns to lowercase
    # the -s option is critical to enable infernal2mapping
    cmd = f"ali-pfam-lowercase-rf-gap-columns.pl -s {temp_stk} > {temp_pfam_stk}"
    result = runner.run(cmd)
    if result:
        raise ValueError(
            f"Failed ali-pfam-lowercase-rf-gap-columns for {seq_id} {model_id}"
        )

    cmd = f"ali-pfam-lowercase-rf-gap-columns.pl -s {temp_stk_original} > {temp_pfam_stk_original}"
    result = runner.run(cmd)
    if result:
        raise ValueError(
            f"Failed ali-pfam-lowercase-rf-gap-columns for {seq_id} {model_id}"
        )

    if not constraint:
        shared.remove_large_insertions_pfam_stk(temp_pfam_stk)
        shared.remove_large_insertions_pfam_stk(temp_pfam_stk_original)

    # convert stockholm to aligned fasta with WUSS secondary structure
    cmd = f"ali-pfam-sindi2dot-bracket.pl -l -n -w -a -c {temp_pfam_stk} > {temp_afa}"
    result = runner.run(cmd)
    if result:
        raise ValueError(f"Failed ali-pfam-sindi2dot-bracket for {seq_id} {model_id}")

    cmd = (
        f"ali-pfam-sindi2dot-bracket.pl -l -n -w -a -c {temp_pfam_stk_original} > "
        f"{temp_afa_original}"
    )
    result = runner.run(cmd)
    if result:
        raise ValueError(f"Failed ali-pfam-sindi2dot-bracket for {seq_id} {model_id}")

    # add original, non-depaired secondary structure
    with open(temp_afa_original) as f_temp_afa_original:
        lines = f_temp_afa_original.readlines()
        ss_cons_original = lines[5]
    with open(temp_afa, "a") as f_temp_afa:
        f_temp_afa.write(">SS_cons_original\n")
        f_temp_afa.write(ss_cons_original)

    # generate traveler infernal mapping file
    infernal_mapping_failed = True
    cmd = f"python3 /rna/traveler/utils/infernal2mapping.py -i {temp_afa} > {temp_map}"
    infernal_mapping_failed = runner.run(cmd)

    if rna_type == "gtrnadb":
        result_base = os.path.join(
            output_folder, seq_id.replace("/", "-") + "-" + domain + "_" + isotype
        )
    else:
        result_base = os.path.join(
            output_folder,
            f"{seq_id.replace('/', '_')}-{model_id}",
        )

    # convert stockholm to fasta with dot bracket secondary structure
    cmd = f"ali-pfam-sindi2dot-bracket.pl {temp_pfam_stk} > {result_base}.fasta"
    result = runner.run(cmd)
    if result:
        rprint(f"Failed esl-pfam-sindi2dot-bracket for {seq_id} {model_id}")
        return

    if constraint:
        shared.fold_insertions(
            f"{result_base}.fasta",
            exclusion,
            rna_type,
            temp_pfam_stk,
            model_id,
            fold_type,
        )
    elif exclusion:
        rprint("Exclusion ignored, enable --constraint to add exclusion file")

    if rna_type == "crw":
        traveler_params = (
            f"--template-structure {template_layout}/{model_id}.ps "
            f"{template_structure}/{model_id}.fasta"
        )
    elif rna_type in ["rfam", "local_data"]:
        traveler_params = (
            f"--template-structure --file-format traveler "
            f"{template_layout} {template_structure} "
        )
    elif rna_type == "gtrnadb":
        traveler_params = (
            f"--template-structure --file-format traveler "
            f"{template_layout} {template_structure} "
            f'--numbering "13,26" -l '
        )
    else:
        traveler_params = (
            f"--template-structure --file-format traveler "
            f"{template_layout}/{model_id}.tr "
            f"{template_structure}/{model_id}.fasta"
        )

    log = result_base + ".log"
    cmd = (
        "traveler --verbose "
        f"--target-structure {result_base}.fasta {traveler_params} "
        f"--draw {temp_map} {result_base} > {log}"
    )
    if not infernal_mapping_failed:
        traveler_failed = runner.run(cmd)

    if infernal_mapping_failed or traveler_failed:
        rprint("Traveler with Infernal mapping failed:")
        rprint(cmd)
        rprint("Repeating using Traveler mapping:")
        cmd = (
            "traveler --verbose "
            f"--target-structure {result_base}.fasta {traveler_params} "
            f"--all {result_base} > {log}"
        )
        runner.run(cmd)

    overlaps = 0
    with open(log) as raw:
        for line in raw:
            match = re.search(r"Overlaps count: (\d+)", line)
            if match:
                if overlaps:
                    rprint("ERROR: Saw too many overlaps")
                    break
                overlaps = int(match.group(1))
    with open(f"{result_base}.overlaps", "w") as out:
        out.write(f"{overlaps}\n")
    if rna_type != "rnasep":
        adjust_font_size(result_base)

    # add metadata to json file
    cmd = (
        f"python3 /rna/traveler/utils/enrich_json.py --input-json {result_base}.colored.json "
        f"--input-data {temp_post_prob} --output {result_base}.enriched.json"
    )
    result = runner.run(cmd)

    # add colors
    if result == 0:
        cmd = (
            f"python3 /rna/traveler/utils/json2svg.py -p /rna/r2dt/utils/colorscheme.json "
            f"-i {result_base}.enriched.json -o {result_base}.enriched.svg"
        )
        runner.run(cmd)

    # clean up
    os.remove(temp_fasta)
    os.remove(temp_sto)
    os.remove(temp_depaired)
    os.remove(temp_stk)
    os.remove(temp_afa)
    os.remove(temp_map)
    if rna_type == "rfam":
        os.remove(temp_sto_unfiltered)
        os.remove(temp_acc_list)


def adjust_font_size(result_base):
    """
    Decrease font-size for large diagrams.
    """
    filenames = [result_base + ".colored.svg", result_base + ".svg"]
    for filename in filenames:
        if not os.path.exists(filename):
            continue
        content = (
            Path(filename).read_text().replace("font-size: 7px;", "font-size: 4px;")
        )
        Path(filename).write_text(content)


# pylint: disable-next=too-many-arguments
def visualise_trna(
    domain, isotype, fasta_input, output_folder, constraint, exclusion, fold_type
):
    """A wrapper for visualising multiple tRNA sequences in a FASTA file."""
    filename = "headers.txt"
    os.makedirs(output_folder, exist_ok=True)

    if not os.path.exists(f"{fasta_input}.ssi"):
        cmd = f"esl-sfetch --index {fasta_input}"
        runner.run(cmd)

    cmd = f"grep '>' {fasta_input} > {filename}"
    runner.run(cmd)

    with open(filename) as f_headers:
        for _, line in enumerate(f_headers):
            seq_id = line.split(" ", 1)[0].replace(">", "").strip()
            rprint(seq_id)
            visualise(
                "gtrnadb",
                fasta_input,
                output_folder,
                seq_id,
                None,
                constraint,
                exclusion,
                fold_type,
                domain,
                isotype,
                None,
                None,
            )
    file_path = Path(filename)
    file_path.unlink(missing_ok=True)
