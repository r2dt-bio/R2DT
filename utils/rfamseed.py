"""A module for managing Rfam seed files."""

import gzip
import os
import re
import shutil
import tempfile
from pathlib import Path

import requests

from . import config
from .runner import runner


class RfamSeed:
    """A class for managing Rfam seed files."""

    def __init__(self) -> None:
        self.seed_archive = self._get_seed_archive()
        self.seed_index = self.seed_archive.with_suffix(".seed.ssi")

    def _get_seed_archive(self) -> Path:
        """Get a path to an Rfam seed alignment archive."""
        return Path(config.CM_LIBRARY) / "rfam" / "Rfam.seed"

    def _index_seed_archive(self) -> None:
        """Index Rfam seed alignment archive."""
        self.seed_index.unlink(missing_ok=True)
        cmd = f"esl-afetch --index {self.seed_archive}"
        runner.run(cmd)

    def download_rfam_seed_archive(self):
        """Get Rfam seed alignment archive
        and store in the downloadable part of the template library."""
        url = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz"
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for HTTP errors
        # Write the file to disk
        with open(f"{self.seed_archive}.gz", "wb") as out_file:
            for chunk in response.iter_content(chunk_size=8192):
                out_file.write(chunk)
        # Decompress the file
        with gzip.open(f"{self.seed_archive}.gz", "rb") as f_in:
            with open(self.seed_archive, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        # Remove the .gz file after decompression
        Path(f"{self.seed_archive}.gz").unlink()
        # Convert to pfam format
        self.convert_rfam_seed_to_pfam()
        # Index the file
        self._index_seed_archive()
        return self

    def convert_rfam_seed_to_pfam(self):
        """Convert to pfam format for ease of parsing
        (each sequence is on 1 line, not split by width)."""
        pfam_seed = Path(f"{self.seed_archive}.pfam")
        cmd = f"esl-reformat pfam {self.seed_archive} > {pfam_seed}"
        runner.run(cmd)
        self.seed_archive.unlink()
        # Rename the pfam file
        pfam_seed.rename(self.seed_archive)

    def get_seed_filename(self, rfam_acc) -> Path:
        """Get a path to an Rfam seed alignment given an accession."""
        filename = Path(tempfile.gettempdir()) / "rfam_seed" / f"{rfam_acc}.seed"
        if not filename.parent.exists():
            filename.parent.mkdir(parents=True, exist_ok=True)
        return filename

    def download_rfam_seed(self, rfam_acc):
        """Fetch Rfam seed alignment using the API."""
        seed_filename = self.get_seed_filename(rfam_acc)
        if seed_filename.exists():
            return seed_filename
        url = f"https://rfam.org/family/{rfam_acc}/alignment"
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for HTTP errors
        seed_filename.write_text(response.text)
        return seed_filename

    def get_rfam_seed(self, rfam_acc) -> Path:
        """Get a path to an Rfam seed alignment given an accession."""
        seed_filename = self.get_seed_filename(rfam_acc)
        if seed_filename.exists():
            return seed_filename
        if self.seed_archive.exists():
            cmd = f"esl-afetch {self.seed_archive} {rfam_acc} > {seed_filename}"
            runner.run(cmd)
        else:
            self.download_rfam_seed(rfam_acc)
        if not seed_filename.exists():
            raise FileNotFoundError(f"Rfam seed alignment not found in {seed_filename}")
        return seed_filename

    def get_no_structure_file(self):
        """Create a file listing all Rfam accessions that have no secondary structure."""
        if not self.seed_archive.exists():
            self.download_rfam_seed_archive()
        # pylint: disable=consider-using-with
        tmpfile = tempfile.NamedTemporaryFile(delete=False)
        cmd = f"grep -E '(#=GF AC|#=GC SS_cons)' {self.seed_archive} > {tmpfile.name}"
        runner.run(cmd)
        no_structure = []
        with open(tmpfile.name, "r") as f_seed:
            for line in f_seed.readlines():
                match = re.match(r"^#=GF\s+AC\s+(RF\d{5})$", line)
                if match:
                    rfam_acc = match.group(1)
                if line.startswith("#=GC SS_cons"):
                    parts = line.split()
                    if "<" not in parts[2] and ">" not in parts[2]:
                        no_structure.append(rfam_acc)
        no_structure_filename = os.path.join(config.RFAM_DATA, "no_structure.txt")
        with open(no_structure_filename, "w") as f_list:
            for rfam_acc in no_structure:
                f_list.write(rfam_acc + "\n")
        tmpfile.close()
        os.unlink(tmpfile.name)
