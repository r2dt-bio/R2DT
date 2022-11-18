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
from . import config
from . import gtrnadb
from . import shared
from . import rfam


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
    elif rna_type.lower() == "gtrnadb":
        model_id = domain + "_" + isotype
    else:
        print("Please specify RNA type")
        return

    temp_fasta = filename_template.replace("type", "fasta")
    temp_sto = filename_template.replace("type", "sto")
    temp_depaired = filename_template.replace("type", "depaired")
    temp_stk = filename_template.replace("type", "stk")
    temp_post_prob = filename_template.replace("type", "post_prob")
    temp_pfam_stk = filename_template.replace("type", "pfam_stk")
    temp_afa = filename_template.replace("type", "afa")
    temp_map = filename_template.replace("type", "map")

    # get sequence from fasta file
    seq_range = f"-c {start}..{end}" if start and end else ""
    cmd = f"esl-sfetch {seq_range} {fasta_input} {seq_id} > {temp_fasta}"
    result = os.system(cmd)
    if result:
        raise ValueError(f"Failed esl-sfetch for: {seq_id}")

    # check that the model exists
    if rna_type == "rfam":
        model_path = rfam.get_rfam_cm(model_id)
        template_layout = rfam.get_traveler_template_xml(model_id)
        template_structure = rfam.get_traveler_fasta(model_id)
        # download seed alignment and list its accessions
        rfam_seed = rfam.download_rfam_seed(model_id)
        cmd = f"esl-alistat --list {temp_acc_list} {rfam_seed} > /dev/null"
        os.system(cmd)
    elif rna_type == "gtrnadb":
        model_path = gtrnadb.get_trnascan_cm(domain, isotype)
        template_layout = gtrnadb.get_traveler_template_xml(domain, isotype)
        template_structure = gtrnadb.get_traveler_fasta(domain, isotype)
    else:
        model_path = os.path.join(cm_library, model_id + ".cm")
        if not os.path.exists(model_path):
            print(f"Model not found {model_path}")
            return

    # align sequence to the model
    cm_options = ["", "--mxsize 2048 --maxtau 0.49"]
    for options in cm_options:
        if rna_type == "rfam":
            cmd = f"cmalign --mapali {rfam_seed} --mapstr {options} {model_path} {temp_fasta} > {temp_sto_unfiltered}"
        else:
            cmd = f"cmalign {options} {model_path} {temp_fasta} > {temp_sto}"
        result = os.system(cmd)
        if not result:
            break
    else:
        print(f"Failed cmalign of {seq_id} to {model_id}")
        return

    if rna_type == "rfam":
        cmd = f"/rna/easel/miniapps/esl-alimanip --seq-r {temp_acc_list} {temp_sto_unfiltered} > {temp_sto}"
        os.system(cmd)

    # remove non-canonical Watson-Crick basepairs (e.g. C:A in URS000008DB9C_7227)
    cmd = f"esl-alidepair.pl --nc 0.5 {temp_sto} {temp_depaired} > /dev/null"
    result = os.system(cmd)
    if result:
        print(f"Failed esl-alidepair for {seq_id}")

    has_conserved_structure = False
    with open(temp_sto, "r", encoding="utf-8") as f_stockholm:
        for line in f_stockholm.readlines():
            if line.startswith("#=GC SS_cons"):
                if "<" in line:
                    has_conserved_structure = True
                else:
                    print("This RNA does not have a conserved structure")
                break
    if not has_conserved_structure:
        return

    # impose consensus secondary structure and convert to pfam format
    cmd = f"/rna/easel/miniapps/esl-alimanip --rna --sindi --outformat pfam {temp_depaired} > {temp_stk}"
    result = os.system(cmd)
    if result:
        print(f"Failed esl-alimanip for {seq_id} {model_id}")
        return

    # store posterior probabilities in tsv file
    shared.get_infernal_posterior_probabilities(temp_stk, temp_post_prob)

    # convert nts that are in RF-gap columns to lowercase
    cmd = f"ali-pfam-lowercase-rf-gap-columns.pl {temp_stk} > {temp_pfam_stk}"
    result = os.system(cmd)
    if result:
        raise ValueError(
            f"Failed ali-pfam-lowercase-rf-gap-columns for {seq_id} {model_id}"
        )

    if not constraint:
        shared.remove_large_insertions_pfam_stk(temp_pfam_stk)

    # convert stockholm to aligned fasta with WUSS secondary structure
    cmd = f"ali-pfam-sindi2dot-bracket.pl -l -n -w -a -c {temp_pfam_stk} > {temp_afa}"
    result = os.system(cmd)
    if result:
        raise ValueError(f"Failed ali-pfam-sindi2dot-bracket for {seq_id} {model_id}")

    # generate traveler infernal mapping file
    infernal_mapping_failed = True
    cmd = f"python3 /rna/traveler/utils/infernal2mapping.py -i {temp_afa} > {temp_map}"
    infernal_mapping_failed = os.system(cmd)

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
    result = os.system(cmd)
    if result:
        print(f"Failed esl-pfam-sindi2dot-bracket for {seq_id} {model_id}")
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
        print("Exclusion ignored, enable --constraint to add exclusion file")

    if rna_type == "crw":
        traveler_params = f"--template-structure {template_layout}/{model_id}.ps {template_structure}/{model_id}.fasta"
    elif rna_type == "rfam":
        traveler_params = f"--template-structure --file-format traveler {template_layout} {template_structure} "
    elif rna_type == "gtrnadb":
        traveler_params = f'--template-structure --file-format traveler {template_layout} {template_structure} --numbering "13,26" -l '
    else:
        traveler_params = f"--template-structure --file-format traveler {template_layout}/{model_id}.tr {template_structure}/{model_id}.fasta"

    log = result_base + ".log"
    cmd = (
        "traveler --verbose "
        f"--target-structure {result_base}.fasta {traveler_params} "
        f"--draw {temp_map} {result_base} > {log}"
    )
    if not infernal_mapping_failed:
        print(cmd)
        traveler_failed = os.system(cmd)

    if infernal_mapping_failed or traveler_failed:
        print("Repeating using Traveler mapping")
        cmd = (
            "traveler --verbose "
            f"--target-structure {result_base}.fasta {traveler_params} "
            f"--all {result_base} > {log}"
        )
        print(cmd)
        os.system(cmd)

    overlaps = 0
    with open(log, "r", encoding="utf-8") as raw:
        for line in raw:
            match = re.search(r"Overlaps count: (\d+)", line)
            if match:
                if overlaps:
                    print("ERROR: Saw too many overlaps")
                    break
                overlaps = int(match.group(1))
    with open(f"{result_base}.overlaps", "w", encoding="utf-8") as out:
        out.write(f"{overlaps}\n")
    if rna_type != "rnasep":
        adjust_font_size(result_base)

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
        cmd = f"""sed -i 's/font-size: 7px;/font-size: 4px;/' {filename}"""
        os.system(cmd)
