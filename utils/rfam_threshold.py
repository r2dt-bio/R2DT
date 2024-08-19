"""
A Python class for finding Rfam model-specific gathering thresholds.
"""

import tempfile
from dataclasses import dataclass
from pathlib import Path

from utils import config


@dataclass
class RfamThreshold:
    """
    A dataclass to store Rfam model-specific gathering thresholds.
    """

    rfam_acc: str
    name: str
    description: str
    gathering_threshold: float
    trusted_cutoff: float
    noise_cutoff: float


class RfamThresholds:
    """
    A class to store Rfam model-specific gathering thresholds.
    """

    def __init__(self):
        self.cm_filename = Path(config.CM_LIBRARY) / "rfam" / "all.cm"
        self.cached_filename = Path(tempfile.gettempdir()) / "rfam_thresholds.tsv"
        self.thresholds = {}
        self.load()

    def load(self):
        """
        Load the Rfam model-specific gathering thresholds.
        """
        if self.cached_filename.exists():
            self._load_cached()
        else:
            self._load_cm()

    def _load_cached(self):
        with open(self.cached_filename) as file:
            for line in file:
                (
                    rfam_acc,
                    name,
                    description,
                    gathering_threshold,
                    trusted_cutoff,
                    noise_cutoff,
                ) = line.strip().split("\t")
                self.thresholds[name] = RfamThreshold(
                    rfam_acc,
                    name,
                    description,
                    float(gathering_threshold),
                    float(trusted_cutoff),
                    float(noise_cutoff),
                )

    def _load_cm(self):
        with open(self.cm_filename) as file:
            for line in file:
                if line.startswith("ACC"):
                    rfam_acc = line.strip().split()[1]
                elif line.startswith("NAME"):
                    name = line.strip().split()[1]
                elif line.startswith("DESC"):
                    description = " ".join(line.strip().split()[1:])
                elif line.startswith("GA"):
                    gathering_threshold = float(line.strip().split()[1])
                elif line.startswith("TC"):
                    trusted_cutoff = float(line.strip().split()[1])
                elif line.startswith("NC"):
                    noise_cutoff = float(line.strip().split()[1])
                    self.thresholds[name] = RfamThreshold(
                        rfam_acc,
                        name,
                        description,
                        gathering_threshold,
                        trusted_cutoff,
                        noise_cutoff,
                    )

        with open(self.cached_filename, "w") as file:
            for _, threshold in self.thresholds.items():
                file.write(
                    "\t".join(
                        [
                            threshold.rfam_acc,
                            threshold.name,
                            threshold.description,
                            str(threshold.gathering_threshold),
                            str(threshold.trusted_cutoff),
                            str(threshold.noise_cutoff),
                        ]
                    )
                    + "\n"
                )

    def get_threshold(self, rfam_name):
        """
        Get the Rfam model-specific gathering thresholds.
        """
        family = self.thresholds.get(rfam_name, None)
        if family is None:
            return None
        return family.gathering_threshold
