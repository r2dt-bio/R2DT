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

import contextlib
import gzip
import re
import shutil
import tempfile
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


# ============================================================================
# Local structure file utilities
# ============================================================================

# Valid structure file extensions (including compressed variants)
STRUCTURE_EXTENSIONS = {".pdb", ".cif", ".pdb.gz", ".cif.gz", ".ent", ".ent.gz"}


def is_gzip_file(file_path: Path) -> bool:
    """
    Check if a file is gzip-compressed by examining magic bytes.

    Args:
        file_path: Path to the file to check.

    Returns:
        True if file starts with gzip magic bytes (0x1f 0x8b).
    """
    try:
        with open(file_path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except (OSError, IOError):
        return False


def get_structure_format(file_path: Path) -> Optional[str]:
    """
    Determine structure format from file extension.

    Handles compressed files by stripping .gz extension first.

    Args:
        file_path: Path to the structure file.

    Returns:
        Format string ("pdb" or "cif") or None if unknown.
    """
    name = file_path.name.lower()

    # Strip .gz if present
    if name.endswith(".gz"):
        name = name[:-3]

    # Check extension
    if name.endswith(".cif"):
        return "cif"
    if name.endswith((".pdb", ".ent")):
        return "pdb"

    return None


# pylint: disable=too-many-return-statements
def validate_structure_file(file_path: Path) -> Tuple[bool, str, Optional[str]]:
    """
    Validate a local PDB/mmCIF structure file.

    Checks that the file exists, is readable, and contains valid structure data.
    Supports gzip-compressed files (.gz).

    Args:
        file_path: Path to the structure file.

    Returns:
        Tuple of (is_valid, format_or_error, error_message).
        On success: (True, "pdb"|"cif", None)
        On failure: (False, "", error_message)
    """
    path = Path(file_path)

    # Check file exists
    if not path.exists():
        return False, "", f"File not found: {path}"

    if not path.is_file():
        return False, "", f"Not a file: {path}"

    # Check file size (catch empty files)
    if path.stat().st_size == 0:
        return False, "", "File is empty"

    # Detect compression
    is_compressed = is_gzip_file(path)

    # Open file appropriately
    try:
        if is_compressed:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
                first_lines = [f.readline() for _ in range(10)]
        else:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                first_lines = [f.readline() for _ in range(10)]
    except gzip.BadGzipFile:
        return False, "", "Corrupted gzip file"
    except OSError as e:
        return False, "", f"Error reading file: {e}"

    first_content = "".join(first_lines)

    # Detect format from content
    detected_format = None

    # CIF files start with 'data_'
    if first_content.lstrip().startswith("data_"):
        detected_format = "cif"
    # PDB files have characteristic record types
    elif any(
        first_content.startswith(kw)
        for kw in ["HEADER", "ATOM", "HETATM", "REMARK", "TITLE", "COMPND", "SOURCE"]
    ):
        detected_format = "pdb"
    else:
        # Try extension-based detection as fallback
        detected_format = get_structure_format(path)

    if not detected_format:
        preview = first_content[:100].replace("\n", " ")
        return False, "", f"Unknown structure format. Content preview: {preview}..."

    return True, detected_format, None


def is_local_structure_file(input_str: str) -> bool:
    """
    Check if input string looks like a local structure file path.

    Args:
        input_str: Input string (could be PDB ID or file path).

    Returns:
        True if input appears to be a file path to a structure file.
    """
    path = Path(input_str)

    # Check if it's an existing file
    if path.exists() and path.is_file():
        # Check extension
        name = path.name.lower()
        for ext in STRUCTURE_EXTENSIONS:
            if name.endswith(ext):
                return True

    # Check if path has structure extension even if file doesn't exist yet
    # (useful for error messages)
    name = path.name.lower()
    for ext in STRUCTURE_EXTENSIONS:
        if name.endswith(ext):
            return True

    return False


@contextlib.contextmanager
def open_structure_file(file_path: Path, mode: str = "r"):
    """
    Context manager to open a structure file, handling gzip compression.

    Automatically detects and decompresses gzip files.

    Args:
        file_path: Path to the structure file.
        mode: File mode ('r' for text, 'rb' for binary).

    Yields:
        File handle (text or binary depending on mode).
    """
    path = Path(file_path)
    is_compressed = is_gzip_file(path)

    if is_compressed:
        gzip_mode = "rt" if mode == "r" else "rb"
        with gzip.open(path, gzip_mode) as f:
            yield f
    else:
        with open(path, mode) as f:
            yield f


class DecompressedStructureFile:
    """
    Context manager that provides a decompressed structure file path.

    For gzip files, decompresses to a temporary file that is automatically
    cleaned up. For uncompressed files, returns the original path.

    This is useful for tools that require a file path on disk (like RNAView)
    rather than a file handle.
    """

    def __init__(self, file_path: Path):
        """
        Initialize with a structure file path.

        Args:
            file_path: Path to the structure file (may be gzip-compressed).
        """
        self.original_path = Path(file_path)
        self.temp_file = None
        self.decompressed_path = None

    def __enter__(self) -> Path:
        """
        Decompress file if needed and return path to readable file.

        Returns:
            Path to the decompressed structure file.
        """
        if not is_gzip_file(self.original_path):
            # Not compressed, return original path
            self.decompressed_path = self.original_path
            return self.decompressed_path

        # Determine output extension
        name = self.original_path.name
        if name.endswith(".gz"):
            decompressed_name = name[:-3]
        else:
            decompressed_name = name

        # Create temporary file with correct extension
        suffix = Path(decompressed_name).suffix or ".pdb"
        self.temp_file = tempfile.NamedTemporaryFile(
            mode="wb", suffix=suffix, delete=False
        )

        # Decompress
        with gzip.open(self.original_path, "rb") as f_in:
            shutil.copyfileobj(f_in, self.temp_file)

        self.temp_file.close()
        self.decompressed_path = Path(self.temp_file.name)
        return self.decompressed_path

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Clean up temporary file if created."""
        if self.temp_file is not None:
            try:
                Path(self.temp_file.name).unlink()
            except OSError:
                pass
        return False
