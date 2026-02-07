"""
Utilities for fetching PDB/mmCIF structure files from RCSB.

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

import re
from pathlib import Path
from typing import Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

# RCSB download URLs
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"

# Request timeout in seconds
DEFAULT_TIMEOUT = 30


def validate_pdb_id(pdb_id: str) -> bool:
    """
    Validate PDB ID format.

    Traditional PDB IDs are 4 characters (e.g., 1S72).
    Extended PDB IDs can be longer alphanumeric strings.

    Args:
        pdb_id: The PDB identifier to validate.

    Returns:
        True if valid, False otherwise.
    """
    # Traditional 4-character PDB ID or extended alphanumeric
    pattern = r"^[A-Za-z0-9]{4,}$"
    return bool(re.match(pattern, pdb_id))


def _fetch_url(url: str, timeout: int = DEFAULT_TIMEOUT) -> Optional[bytes]:
    """
    Fetch content from a URL using urllib.

    Args:
        url: The URL to fetch.
        timeout: Request timeout in seconds.

    Returns:
        Content as bytes, or None if not available.
    """
    try:
        request = Request(url, headers={"User-Agent": "R2DT/1.0"})
        with urlopen(request, timeout=timeout) as response:
            return response.read()
    except (URLError, HTTPError):
        return None


def fetch_pdb(pdb_id: str, timeout: int = DEFAULT_TIMEOUT) -> Optional[bytes]:
    """
    Fetch PDB format file from RCSB.

    Args:
        pdb_id: The PDB identifier (e.g., "1S72").
        timeout: Request timeout in seconds.

    Returns:
        File content as bytes, or None if not available.
    """
    url = RCSB_PDB_URL.format(pdb_id=pdb_id.upper())
    return _fetch_url(url, timeout)


def fetch_cif(pdb_id: str, timeout: int = DEFAULT_TIMEOUT) -> Optional[bytes]:
    """
    Fetch mmCIF format file from RCSB.

    Args:
        pdb_id: The PDB identifier (e.g., "1S72").
        timeout: Request timeout in seconds.

    Returns:
        File content as bytes, or None if not available.
    """
    url = RCSB_CIF_URL.format(pdb_id=pdb_id.upper())
    return _fetch_url(url, timeout)


def download_structure(
    pdb_id: str,
    output_dir: Path,
    prefer_format: str = "pdb",
    timeout: int = DEFAULT_TIMEOUT,
) -> Tuple[Optional[Path], Optional[str]]:
    """
    Download structure from RCSB, with automatic fallback.

    Tries to download in the preferred format first. If unavailable,
    falls back to the alternative format.

    Args:
        pdb_id: The PDB identifier (e.g., "1S72" or "9FN3").
        output_dir: Directory to save the downloaded file.
        prefer_format: Preferred format ("pdb" or "cif"). Default is "pdb".
        timeout: Request timeout in seconds.

    Returns:
        Tuple of (file_path, format) where format is "pdb" or "cif".
        Returns (None, None) if download fails.

    Raises:
        ValueError: If pdb_id is invalid or prefer_format is not recognized.
    """
    pdb_id = pdb_id.upper()

    if not validate_pdb_id(pdb_id):
        raise ValueError(f"Invalid PDB ID: {pdb_id}")

    if prefer_format not in ("pdb", "cif"):
        raise ValueError(f"Invalid format: {prefer_format}. Must be 'pdb' or 'cif'.")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine fetch order based on preference
    if prefer_format == "pdb":
        fetch_order = [("pdb", fetch_pdb), ("cif", fetch_cif)]
    else:
        fetch_order = [("cif", fetch_cif), ("pdb", fetch_pdb)]

    for fmt, fetch_func in fetch_order:
        content = fetch_func(pdb_id, timeout)
        if content:
            file_path = output_dir / f"{pdb_id}.{fmt}"
            file_path.write_bytes(content)
            return file_path, fmt

    return None, None
