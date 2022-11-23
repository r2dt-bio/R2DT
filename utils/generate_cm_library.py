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


def convert_bpseq_to_fasta(bpseq):
    """Use a Traveler script to convert from BPSEQ to FASTA."""
    fasta = bpseq.replace(".bpseq", ".fasta")
    if not os.path.exists(fasta):
        cmd = f"python /rna/traveler/utils/bpseq2fasta.py -i {bpseq} -o {fasta}"
        os.system(cmd)
    return fasta


def break_pseudoknots(fasta):
    """Remove pseudoknots using RNAStructure."""
    fasta_no_knots = fasta.replace("-with-knots.fasta", ".fasta")
    if not os.path.exists(fasta_no_knots):
        cmd = f"RemovePseudoknots -b {fasta} {fasta_no_knots}"
        os.system(cmd)
    return fasta_no_knots


def convert_fasta_to_stockholm(fasta):
    """Convert fasta to stockholm."""
    stockholm = fasta.replace(".fasta", ".sto")
    model_id = os.path.basename(stockholm).replace(".sto", "")
    if not os.path.exists(stockholm):
        with open(fasta, "r", encoding="utf-8") as f_input:
            with open(stockholm, "w", encoding="utf-8") as f_output:
                lines = f_input.readlines()
                f_output.write("# STOCKHOLM 1.0\n")
                f_output.write("\n")
                f_output.write(f"{model_id.ljust(60)}{lines[1].strip()}\n")
                f_output.write(f"{'#=GC SS_cons'.ljust(60)}{lines[2].strip()}\n")
                f_output.write("//\n")
    return stockholm


def copy_cm_evalues(cm_filename):
    """
    Update covariance files generated from CRW covariance models
    by copying E-values from Rfam CMs.
    """
    rfam_acc = "RF00177"
    example_cm_filename = os.path.join("temp", f"{rfam_acc}.cm")
    if not os.path.exists(example_cm_filename):
        cmd = f"wget -O {rfam_acc}.cm http://rfam.org/family/{rfam_acc}/cm"
        os.system(cmd)
    cmd = (
        f"perl /rna/jiffy-infernal-hmmer-scripts/cm-copy-evalue-parameters.pl "
        f"{rfam_acc}.cm {cm_filename}"
    )
    os.system(cmd)
    os.system(f"rm {cm_filename}.old")


def build_cm(stockholm, cm_library):
    """Build an Infernal covariance model."""
    cm_filename = os.path.join(
        cm_library, os.path.basename(stockholm).replace(".sto", ".cm")
    )
    if not os.path.exists(cm_filename):
        cmd = f"cmbuild {cm_filename} {stockholm}"
        os.system(cmd)
        copy_cm_evalues(cm_filename)
    else:
        print(f"CM already exists {cm_filename}")
    return cm_filename
