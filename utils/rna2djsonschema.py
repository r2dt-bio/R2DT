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
from pathlib import Path

from utils import config
from utils import generate_cm_library as gcl

TEMPLATE_FOLDER = "new"


# pylint: disable=too-many-instance-attributes
class SchemaToTemplate:
    """Convert an RNA 2D JSON Schema file to a Traveler template."""

    def __init__(self, json_file: str):
        self.json_file = Path(json_file)
        self.rna_name = self.json_file.stem
        self.data = self.get_data()
        self.path = self.get_file_location()
        self.sequence = self.get_sequence()
        self.dot_bracket = self.get_dot_bracket()
        self.prime_label_present = False
        self.result = {
            "xml": self.path / f"{self.rna_name}.xml",
            "fasta": self.path / f"{self.rna_name}.fasta",
            "cm": self.path / f"{self.rna_name}.cm",
        }
        self.create_template()

    def __repr__(self) -> str:
        """Return a string representation of the object."""
        return (
            f"SchemaToTemplate("
            f"path={self.path}, "
            f"cm={self.result['cm'].name}, "
            f"fasta={self.result['fasta'].name}, "
            f"xml={self.result['xml'].name}"
            f")"
        )

    def get_data(self) -> dict:
        """Parse RNA 2D JSON Schema file and return the data and other info."""
        with open(self.json_file) as f_json:
            data = json.load(f_json)
        return data

    def get_file_location(self) -> Path:
        """Get the location of the output files."""
        destination = Path(config.DATA) / TEMPLATE_FOLDER / self.rna_name
        destination.mkdir(parents=True, exist_ok=True)
        return destination

    def generate_traveler_xml(self) -> None:
        """Generate a Traveler XML file from an RNA 2D JSON Schema file."""
        molecule = self.data["rnaComplexes"][0]["rnaMolecules"][0]
        with open(self.result["xml"], "w") as f_xml:
            f_xml.write("<structure>\n")
            for nucleotide in molecule["sequence"]:
                if nucleotide["residueName"] in ["5'", "3'"]:
                    self.prime_label_present = True
                    continue
                base = nucleotide["residueName"]
                x_coord = nucleotide["x"]
                y_coord = nucleotide["y"]
                f_xml.write(
                    f'<point x="{x_coord:.2f}" y="{y_coord:.2f}" b="{base}"/>\n'
                )
            f_xml.write("</structure>\n")

    def get_sequence(self) -> str:
        """Get the sequence of the RNA molecule."""
        molecule = self.data["rnaComplexes"][0]["rnaMolecules"][0]
        sequence = []
        for nucleotide in molecule["sequence"]:
            if nucleotide["residueName"] in ["5'", "3'"]:
                self.prime_label_present = True
                continue
            sequence.append(nucleotide["residueName"])
        return "".join(sequence)

    def get_dot_bracket(self) -> str:
        """Get the dot bracket notation of the RNA molecule."""
        dot_bracket = ["."] * len(self.sequence)
        molecule = self.data["rnaComplexes"][0]["rnaMolecules"][0]
        for basepair in molecule["basePairs"]:
            if basepair["basePairType"] != "canonical":
                continue
            if self.prime_label_present:
                basepair["residueIndex1"] -= 1
                basepair["residueIndex2"] -= 1
            dot_bracket[basepair["residueIndex1"]] = "("
            dot_bracket[basepair["residueIndex2"]] = ")"
        return "".join(dot_bracket)

    def generate_traveler_fasta(self) -> None:
        """Generate a Traveler FASTA file from an RNA 2D JSON Schema file."""
        with open(self.result["fasta"], "w", encoding="utf-8") as f_fasta:
            lines = [
                f">{self.rna_name}",
                self.sequence,
                self.dot_bracket,
            ]
            f_fasta.write("\n".join(lines) + "\n")

    def generate_cm(self) -> None:
        """Generate a covariance model from an RNA 2D JSON Schema file."""
        gcl.build_cm(
            gcl.convert_fasta_to_stockholm(str(self.result["fasta"])), self.path
        )

    def create_template(self) -> None:
        """Create a template."""
        self.generate_traveler_xml()
        self.generate_traveler_fasta()
        self.generate_cm()
