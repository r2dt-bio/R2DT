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
import glob


# CRW_FASTA_NO_PSEUDOKNOTS = '/rna/r2dt/data/crw-fasta-no-pseudoknots'


def convert_bpseq_to_fasta(bpseq):
    fasta = bpseq.replace(".bpseq", ".fasta")
    if not os.path.exists(fasta):
        cmd = f"python /rna/traveler/utils/bpseq2fasta.py -i {bpseq} -o {fasta}"
        os.system(cmd)
    return fasta


def break_pseudoknots(fasta):
    fasta_no_knots = fasta.replace("-with-knots.fasta", ".fasta")
    if not os.path.exists(fasta_no_knots):
        cmd = f"RemovePseudoknots -b {fasta} {fasta_no_knots}"
        os.system(cmd)
    return fasta_no_knots


def convert_fasta_to_stockholm(fasta):
    stockholm = fasta.replace(".fasta", ".sto")
    model_id = os.path.basename(stockholm).replace(".sto", "")
    if not os.path.exists(stockholm):
        with open(fasta, "r") as f_input:
            with open(stockholm, "w") as f_output:
                lines = f_input.readlines()
                f_output.write("# STOCKHOLM 1.0\n")
                f_output.write("\n")
                f_output.write(f"{model_id.ljust(60)}{lines[1].strip()}\n")
                f_output.write(
                    f"{'#=GC SS_cons'.ljust(60)}{lines[2].strip()}\n"
                )
                f_output.write("//\n")
    return stockholm


def copy_cm_evalues(cm):
    """
    Update covariance files genenrated from CRW covariance models
    by copying E-values from Rfam CMs.
    """
    if not os.path.exists("RF00177.cm"):
        cmd = "wget -O RF00177.cm http://rfam.org/family/RF00177/cm"
        os.system(cmd)
    cmd = "perl /rna/jiffy-infernal-hmmer-scripts/cm-copy-evalue-parameters.pl RF00177.cm {cm}".format(
        cm=cm
    )
    os.system(cmd)
    os.system(f"rm {cm}.old")


def build_cm(stockholm, cm_library):
    cm = os.path.join(cm_library, os.path.basename(stockholm).replace(".sto", ".cm"))
    if not os.path.exists(cm):
        cmd = f"cmbuild {cm} {stockholm}"
        os.system(cmd)
        copy_cm_evalues(cm)
    else:
        print(f"CM already exists {cm}")
    return cm
