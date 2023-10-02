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

import json
import os


def parse_json_file(json_file):
    """Parse RNA 2D JSON Schema file and return the data and other info."""
    with open(json_file) as f_json:
        data = json.load(f_json)
    rna_name = os.path.basename(json_file).replace(".json", "")
    destination = os.path.join("data", "new", rna_name)
    os.makedirs(destination, exist_ok=True)
    return data, destination, rna_name


def generate_traveler_xml(data, destination, rna_name):
    """Generate a Traveler XML file from an RNA 2D JSON Schema file."""
    xml_template = os.path.join(destination, f"{rna_name}.xml")
    with open(xml_template, "w") as f_xml:
        f_xml.write("<structure>\n")
        for nucleotide in data["rnaComplexes"][0]["rnaMolecules"][0]["sequence"]:
            if nucleotide["residueName"] in ["5'", "3'"]:
                continue
            base = nucleotide["residueName"]
            x_coord = nucleotide["x"]
            y_coord = nucleotide["y"]
            f_xml.write(f'<point x="{x_coord:.2f}" y="{y_coord:.2f}" b="{base}"/>\n')
        f_xml.write("</structure>\n")
    return xml_template


def generate_traveler_fasta(data, destination, rna_name):
    """Generate a Traveler FASTA file from an RNA 2D JSON Schema file."""
    fasta_file = os.path.join(destination, f"{rna_name}.fasta")
    sequence = []
    for nucleotide in data["rnaComplexes"][0]["rnaMolecules"][0]["sequence"]:
        if nucleotide["residueName"] in ["5'", "3'"]:
            continue
        sequence.append(nucleotide["residueName"])
    dot_bracket = ["."] * len(sequence)
    for basepair in data["rnaComplexes"][0]["rnaMolecules"][0]["basePairs"]:
        if basepair["basePairType"] != "canonical":
            continue
        dot_bracket[basepair["residueIndex1"]] = "("
        dot_bracket[basepair["residueIndex2"]] = ")"

    with open(fasta_file, "w", encoding="utf-8") as f_fasta:
        f_fasta.write(f">{rna_name}\n")
        f_fasta.write(f"{''.join(sequence)}\n")
        f_fasta.write(f"{''.join(dot_bracket)}\n")
    return fasta_file
