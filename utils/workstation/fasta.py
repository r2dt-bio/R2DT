"""FASTA helpers for workstation draw jobs."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import List, Optional


@dataclass
class FastaRecord:
    """One FASTA entry (optional third-line structure for templatefree)."""

    seq_id: str
    sequence: str
    structure: Optional[str] = None

    @property
    def length(self) -> int:
        """Number of nucleotides in the sequence."""
        return len(self.sequence)

    def to_fasta_text(self) -> str:
        """Serialise as 2- or 3-line FASTA."""
        lines = [f">{self.seq_id}", self.sequence]
        if self.structure:
            lines.append(self.structure)
        return "\n".join(lines) + "\n"


_ID_SAFE = re.compile(r"[^A-Za-z0-9._-]+")
_STRUCT_RE = re.compile(r"^[.()<>{}[\]A-Za-z]+$")


def safe_seq_id(raw: str, fallback: str = "seq") -> str:
    """Filesystem-safe sequence id (no leading '>')."""
    text = (raw or "").strip()
    if text.startswith(">"):
        text = text[1:].strip()
    if text:
        text = text.split()[0]
    text = _ID_SAFE.sub("_", text).strip("._-")
    return text or fallback


def _looks_like_structure(line: str, seq_len: int) -> bool:
    if len(line) != seq_len:
        return False
    if not _STRUCT_RE.match(line):
        return False
    # Prefer lines that contain secondary-structure brackets.
    return any(ch in line for ch in ".()<>{}[]")


def parse_fasta_records(  # pylint: disable=too-many-branches
    text: str,
) -> List[FastaRecord]:
    """Parse multi-FASTA; support optional third-line dot-bracket per record."""
    if not text or not text.strip():
        raise ValueError("FASTA text is empty")

    blocks: List[tuple] = []
    current_id: Optional[str] = None
    body: List[str] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                blocks.append((current_id, body))
            current_id = safe_seq_id(line, fallback=f"seq{len(blocks) + 1}")
            body = []
            continue
        if current_id is None:
            raise ValueError("FASTA must start with a >header line")
        body.append(line)

    if current_id is not None:
        blocks.append((current_id, body))
    if not blocks:
        raise ValueError("No FASTA records found")

    records: List[FastaRecord] = []
    for seq_id, lines in blocks:
        if not lines:
            raise ValueError(f"Empty sequence for {seq_id}")
        structure = None
        if len(lines) >= 2:
            seq_joined = "".join(lines[:-1])
            seq_joined = re.sub(r"[^A-Za-z]", "", seq_joined).upper()
            last = lines[-1]
            if seq_joined and _looks_like_structure(last, len(seq_joined)):
                sequence = seq_joined
                structure = last
            else:
                sequence = re.sub(r"[^A-Za-z]", "", "".join(lines)).upper()
        else:
            sequence = re.sub(r"[^A-Za-z]", "", lines[0]).upper()
        if not sequence:
            raise ValueError(f"Empty sequence for {seq_id}")
        if structure is not None and len(structure) != len(sequence):
            raise ValueError(
                f"Structure length {len(structure)} != sequence length "
                f"{len(sequence)} for {seq_id}"
            )
        records.append(
            FastaRecord(seq_id=seq_id, sequence=sequence, structure=structure)
        )
    return records
