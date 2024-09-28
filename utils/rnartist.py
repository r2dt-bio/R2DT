"""This module contains a class for managing RNArtist layouts."""

import tempfile
from pathlib import Path

from . import config
from .runner import runner


class RnaArtist:
    """A class for managing RNArtist layouts."""

    def __init__(self, rfam_acc="rnartist", destination=None) -> None:
        """Create RnaArtist object."""
        self.rfam_acc = rfam_acc
        if destination:
            self.destination = Path(destination)
        else:
            self.destination = Path(config.RFAM_DATA) / self.rfam_acc
        if self.rfam_acc == "rnartist":
            self.fasta_file = ""
        else:
            self.fasta_file = self.load_rfam_fasta()
        self.stockholm = self.destination / f"{self.rfam_acc}.sto"
        self.rnartist_xml = self.destination / "rnartist-template.xml"
        if not self.destination.exists():
            self.destination.mkdir(parents=True, exist_ok=True)
        self.seq_label = "rnartist"
        self.rnartist_executable = "/usr/local/bin/rnartist.jar"

    # pylint: disable=import-outside-toplevel
    def load_rfam_fasta(self):
        """Import get_traveler_fasta method from Rfam
        avoiding circular imports."""
        from .rfam import get_traveler_fasta

        return get_traveler_fasta(self.rfam_acc)

    def run(self, rerun=False, detail_level=4) -> None:
        """Create a Traveler XML template based on an
        RNArtist consensus 2D layout."""
        if self.rnartist_xml.exists() and not rerun:
            return
        self.make_stockholm_for_rnartist()
        self.run_rnartist(detail_level)
        self.fix_rnartist_nt_labels()
        self.clean_up()

    def make_stockholm_for_rnartist(self):
        """Take Traveler fasta file and create a Stockholm file for RNArtist."""
        with open(self.fasta_file, "r", encoding="utf-8") as f_in:
            lines = f_in.readlines()
            with open(self.stockholm, "w") as f_sto:
                f_sto.write("# STOCKHOLM 1.0\n")
                f_sto.write(f"{self.seq_label}  {lines[1]}")
                f_sto.write(f"#=GC SS_cons      {lines[2]}")
                f_sto.write("//\n")

    def run_rnartist(self, detail_level):
        """
        Run RNArtist to get a Traveler layout.
        Note that curly braces in the template are escaped with another curly brace.
        """
        template = """
        rnartist {{
            svg {{
                path = "{destination}"
            }}
            traveler {{
                path = "{destination}"
                width = 500.0
                height = 500.0
            }}
            ss {{
                stockholm {{
                    file = "{stockholm}"
                    name = "{seq_label}"
                }}
            }}
            theme {{
                details {{
                    value = {detail_level}
                }}
            }}
        }}
        """.format(
            stockholm=self.stockholm.resolve(),
            destination=self.destination.resolve(),
            seq_label=self.seq_label,
            detail_level=detail_level,
        )
        with tempfile.NamedTemporaryFile(suffix=".kts", delete=False) as temp:
            kts_file = Path(temp.name)
            temp.write(template.encode())
        cmd = f"java -jar {self.rnartist_executable} {kts_file}"
        runner.run(cmd)
        kts_file.unlink()

    def fix_rnartist_nt_labels(self):
        """Replace "X" nucleotides with actual nucleotdes in RNArtist Traveler XML."""
        temp_xml = self.destination / f"{self.rfam_acc}_{self.seq_label}.traveler"
        if not temp_xml.exists():
            print(f"Error: RNartist XML file not found for {self.rfam_acc}")
            return
        with open(self.fasta_file, "r", encoding="utf-8") as f_fasta:
            lines = f_fasta.readlines()
            sequence = lines[1].strip()
        updated_xml = []
        with open(temp_xml, "r", encoding="utf-8") as f_in:
            for index, line in enumerate(f_in.readlines()):
                if line.startswith("<point") and 'b="X"' in line:
                    updated_xml.append(line.replace("X", sequence[index - 1]))
                else:
                    updated_xml.append(line)
        with open(self.rnartist_xml, "w", encoding="utf-8") as f_out:
            f_out.write("".join(updated_xml))
        temp_xml.unlink()

    def clean_up(self):
        """Remove temporary files."""
        self.stockholm.unlink()
        seq_kts = self.destination / f"{self.rfam_acc}_{self.seq_label}.kts"
        seq_kts.unlink(missing_ok=True)
