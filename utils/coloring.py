"""
Assign per-panel colors for stitched Stockholm diagrams.

Three modes are supported:

* **structure** – each ``structure_id`` gets a deterministic colour via
  :pypi:`colorhash`.
* **region** – all structures that share a ``region_id`` get the same
  colour; structures without a region fall back to their own name.
* **config** – colours are read from a user-supplied TSV file.

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

import csv
from pathlib import Path
from typing import Optional

from colorhash import ColorHash  # pylint: disable=import-error


def _load_color_config(config_path: Path) -> tuple[dict[str, str], Optional[str]]:
    """Load a TSV colour configuration file.

    The file must have two tab-separated columns (no header row)::

        5'UTR           steelblue
        core_protein    #457b9d
        IRES            red
        *               gray

    The special key ``*`` sets the default colour for names not listed
    explicitly.  Any value that SVG understands is accepted: named
    colours, hex codes, ``rgb()`` expressions, etc.

    Returns:
        A tuple of (mapping, default_color).  *default_color* is ``None``
        when no ``*`` row is present.
    """
    mapping: dict[str, str] = {}
    default_color: Optional[str] = None

    with open(config_path, encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            # Skip blank lines and comments
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 2:
                continue
            key = row[0].strip()
            color = row[1].strip()
            if key == "*":
                default_color = color
            else:
                mapping[key] = color

    return mapping, default_color


def build_color_map(
    regions: list[dict],
    mode: str,
    config_path: Optional[Path] = None,
) -> list[str]:
    """Return a list of SVG colour strings, one per panel.

    Args:
        regions: The ``processed_regions`` list produced by
            :func:`utils.stockholm.process_stockholm_alignment`.
            Each dict must have ``"name"`` and may have ``"region"``.
        mode: One of ``"structure"``, ``"region"``, or ``"config"``.
        config_path: Path to a TSV colour-config file (required when
            *mode* is ``"config"``).

    Returns:
        A list of SVG colour strings parallel to *regions*.

    Raises:
        ValueError: If *mode* is ``"config"`` and a region name cannot
            be resolved.
    """
    colors: list[str] = []

    if mode == "structure":
        for region in regions:
            colors.append(ColorHash(region["name"]).hex)

    elif mode == "region":
        for region in regions:
            key = region.get("region") or region["name"]
            colors.append(ColorHash(key).hex)

    elif mode == "config":
        if config_path is None:
            raise ValueError("config_path is required when mode is 'config'")
        mapping, default_color = _load_color_config(Path(config_path))

        for region in regions:
            name = region["name"]
            parent = region.get("region") or ""

            if name in mapping:
                colors.append(mapping[name])
            elif parent in mapping:
                colors.append(mapping[parent])
            elif default_color is not None:
                colors.append(default_color)
            else:
                raise ValueError(
                    f"No colour defined for structure '{name}' "
                    f"(region '{parent}') and no default ('*') set "
                    f"in {config_path}"
                )
    else:
        raise ValueError(f"Unknown coloring mode: {mode!r}")

    return colors
