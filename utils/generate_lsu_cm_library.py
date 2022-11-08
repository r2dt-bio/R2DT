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


from generate_cm_library import build_cm, convert_fasta_to_stockholm

here = os.path.realpath(os.path.dirname(__file__))
project_folder = os.path.dirname(here)
ribovision_folder = os.path.join(project_folder, "data", "ribovision-lsu")

BPSEQ_LOCATION = os.path.join(ribovision_folder, "bpseq")
CM_LIBRARY = os.path.join(ribovision_folder, "cms")


def convert_bpseq_to_fasta(bpseq):
    fasta = bpseq.replace(".bpseq", ".fasta")
    if not os.path.exists(fasta):
        cmd = "python /rna/traveler/utils/bpseq2fasta.py -i {bpseq} -o {fasta}".format(
            bpseq=bpseq, fasta=fasta
        )
        os.system(cmd)
    return fasta


def main():
    for bpseq in glob.glob("%s/*.bpseq" % BPSEQ_LOCATION):
        print(os.path.basename(bpseq).replace(".bpseq", ""))
        fasta = convert_bpseq_to_fasta(bpseq)
        stockholm = convert_fasta_to_stockholm(fasta)
        build_cm(stockholm, CM_LIBRARY)
    print("Done")


if __name__ == "__main__":
    main()
