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
import re
from pathlib import Path

from . import config
from .runner import runner

SCORE_CUTOFF = 25
TRNASCAN_MODELS = "/usr/lib/tRNAscan-SE/models"
TRNASCAN_CONF = "/usr/bin/tRNAscan-SE.conf"


def setup():
    """Extract tRNAScan covariance models as separate files."""
    base = os.path.join("usr", "lib", "tRNAscan-SE", "models")
    cm_dbs = {
        "TRNAinf-arch-iso": "A",
        "TRNAinf-bact-iso": "B",
        "TRNAinf-euk-iso": "E",
        "TRNAinf-mito-vert": "M",
    }
    for cm_file, domain in cm_dbs.items():
        path = Path(os.path.join(base, cm_file))
        with path.open("r", encoding="utf-8") as raw:
            for line in raw:
                line = line.strip()
                if line.startswith("NAME"):
                    _, name = re.split(r"\s+", line, maxsplit=1)
                    if "mito" not in str(path):
                        _, isotype = name.split("-", 1)
                    else:
                        isotype = name
                    get_trnascan_cm(domain, isotype)


def verify_anticodon(isotype, anticodon, start, end):
    """
    When multiple models are possible, select the most likely one.

    Currently only Cys, Leu, and Ser TRNAinf-mito-vert have more than one CM/template.

    Leu
    LeuTAA and LeuTAG do not have any structural or length differences.
    Either template can be used.

    Cys
    There is only one anticodon for Cys
    <= 61 bp without D-arm.
    >= 65 with D-arm (typical Cys model).
    in between can be either.

    Ser
    >= 70 bp are usually SerTGA (with a D-arm).
    < 65 are usually SerGCT (without a D-arm).
    in between can be either.
    """
    if anticodon != "NNN":
        return anticodon
    seq_length = max(start, end) - min(start, end) + 1
    if isotype == "Leu":
        anticodon = "TAA" if seq_length % 2 == 0 else "TAG"
    elif isotype == "Ser":
        if seq_length >= 70:
            anticodon = "TGA"
        elif seq_length < 65:
            anticodon = "GCT"
        else:
            anticodon = "GCT" if seq_length % 2 == 0 else "TGA"
    return anticodon


def parse_trnascan_output(filename):
    """
    Sequence           		     tRNA	Bounds	tRNA	Anti	Intron  Bounds	Inf
    Name               	tRNA #	Begin	End	    Type	Codon	Begin	   End	Score	Note
    --------           	------	-----	------	----	-----	-----	----	------	------
    URS0000023412_9606 	1	1 	73	Thr	TGT	0	0	60.2
    """
    data = {}
    with open(filename, "r", encoding="utf-8") as f_trnascan:
        for i, line in enumerate(f_trnascan):
            if i in [0, 1, 2]:
                continue  # skip 3 header lines
            parts = [x.strip() for x in line.split("\t")]
            seq_id, _, start, end, isotype, anticodon, _, _, score, note = parts
            score = float(score)
            start = int(start)
            end = int(end)
            if score < SCORE_CUTOFF:
                continue
            data[seq_id] = {
                "score": score,
                "isotype": isotype,
                "anticodon": verify_anticodon(isotype, anticodon, start, end),
                "note": note.lower(),
                "start": start,
                "end": end,
            }
    return data


def run_trnascan(fasta_input, output_folder, domain):
    """Launch tRNAScan-SE and return parsed results."""
    _, extension = os.path.splitext(fasta_input)
    output_file = os.path.join(
        output_folder,
        domain + "-" + os.path.basename(fasta_input).replace(extension, ".txt"),
    )
    if domain == "M":
        domain = "M vert"
    if not os.path.exists(output_file):
        runner.run(
            f"tRNAscan-SE -c {TRNASCAN_CONF} -q -{domain} -o {output_file} {fasta_input}"
        )
    return parse_trnascan_output(output_file)


def skip_trna(entry):
    """
    Some tRNAs should not be drawn and need to be skipped.
    """
    return "pseudo" in entry["note"] or entry["isotype"] in ["Undet", "Sup"]


def classify_trna_sequences(fasta_input, output_folder):
    """Run tRNAScan-SE 2.0 and select the matching model."""
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    mito_vert = run_trnascan(fasta_input, output_folder, "M")
    bacteria = run_trnascan(fasta_input, output_folder, "B")
    archaea = run_trnascan(fasta_input, output_folder, "A")
    eukaryotes = run_trnascan(fasta_input, output_folder, "E")

    rna_ids = set()
    rna_ids.update(
        list(mito_vert.keys()),
        list(bacteria.keys()),
        list(archaea.keys()),
        list(eukaryotes.keys()),
    )

    data = []
    for rna_id in rna_ids:
        values = [
            bacteria[rna_id]["score"] if rna_id in bacteria else -1000,
            archaea[rna_id]["score"] if rna_id in archaea else -1000,
            eukaryotes[rna_id]["score"] if rna_id in eukaryotes else -1000,
            mito_vert[rna_id]["score"] if rna_id in mito_vert else -1000,
        ]
        maximum = values.index(max(values))
        if maximum == 0:
            if skip_trna(bacteria[rna_id]):
                continue
            bacteria[rna_id]["domain"] = "B"
            bacteria[rna_id]["id"] = rna_id
            data.append(bacteria[rna_id])
        elif maximum == 1:
            if skip_trna(archaea[rna_id]):
                continue
            archaea[rna_id]["domain"] = "A"
            archaea[rna_id]["id"] = rna_id
            data.append(archaea[rna_id])
        elif maximum == 2:
            if skip_trna(eukaryotes[rna_id]):
                continue
            eukaryotes[rna_id]["domain"] = "E"
            eukaryotes[rna_id]["id"] = rna_id
            data.append(eukaryotes[rna_id])
        elif maximum == 3:
            mito_vert[rna_id]["domain"] = "M"
            mito_vert[rna_id]["id"] = rna_id
            if mito_vert[rna_id]["isotype"] in ["Leu", "Ser"]:
                # LeuTAA, LeuTAG, SerGCT, SerTGA
                mito_vert[rna_id]["isotype"] = (
                    mito_vert[rna_id]["isotype"] + mito_vert[rna_id]["anticodon"]
                )
            data.append(mito_vert[rna_id])

    with open(os.path.join(output_folder, "hits.txt"), "w", encoding="utf-8") as f_out:
        for entry in data:
            f_out.write(f"{entry['id']}\t{entry['domain']}_{entry['isotype']}\tPASS\n")
    return data


def get_trnascan_cm(domain, isotype):
    """
    Fetch a domain-specific isotype covariance model as a separate file.
    """
    if not os.path.exists(config.GTRNADB_CM_LIBRARY):
        os.mkdir(config.GTRNADB_CM_LIBRARY)
    cm_output = Path(config.GTRNADB_CM_LIBRARY) / f"{domain}_{isotype}.cm"
    if cm_output.exists():
        return str(cm_output)

    cm_library = Path(TRNASCAN_MODELS)
    if domain == "A":
        cm_library = cm_library / "TRNAinf-arch-iso"
        cm_name = "arch-" + isotype
    elif domain == "B":
        cm_library = cm_library / "TRNAinf-bact-iso"
        cm_name = "bact-" + isotype
    elif domain == "E":
        cm_library = cm_library / "TRNAinf-euk-iso"
        cm_name = "euk-" + isotype
    elif domain == "M":
        cm_library = cm_library / "TRNAinf-mito-vert"
        cm_name = isotype
    else:
        raise ValueError(f"Unknown domain: {domain}")

    result = runner.run(f"cmfetch -o {cm_output} {cm_library} {cm_name}")
    if result:
        if cm_output.exists():
            cm_output.unlink()
        cm_output = None
    return cm_output


def get_traveler_template_xml(domain, isotype):
    """Get Traveler template with coordinates."""
    if domain == "A":
        return os.path.join(
            config.GTRNADB_ARCH, f"arch-{isotype}-traveler-template.xml"
        )
    if domain == "B":
        return os.path.join(
            config.GTRNADB_BACT, f"bact-{isotype}-traveler-template.xml"
        )
    if domain == "M":
        if "Leu" in isotype or "Ser" in isotype:
            isotype = isotype[0:3] + "_" + isotype[3:6]
        return os.path.join(
            config.GTRNADB_MITO, f"mito_vert_{isotype}-traveler-template.xml"
        )
    if domain == "E":
        return os.path.join(config.GTRNADB_EUK, f"euk-{isotype}-traveler-template.xml")
    raise ValueError(f"Unknown domain {domain}")


def get_traveler_fasta(domain, isotype):
    """Get Traveler structure file."""
    if domain == "A":
        return os.path.join(config.GTRNADB_ARCH, f"arch-{isotype}-traveler.fasta")
    if domain == "B":
        return os.path.join(config.GTRNADB_BACT, f"bact-{isotype}-traveler.fasta")
    if domain == "M":
        if "Leu" in isotype or "Ser" in isotype:
            isotype = isotype[0:3] + "_" + isotype[3:6]
        return os.path.join(config.GTRNADB_MITO, f"mito_vert_{isotype}-traveler.fasta")
    if domain == "E":
        return os.path.join(config.GTRNADB_EUK, f"euk-{isotype}-traveler.fasta")
    raise ValueError(f"Unknown domain {domain}")
