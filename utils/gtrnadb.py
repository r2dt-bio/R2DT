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
import subprocess as sp
from pathlib import Path

from . import config
from . import ribovision


def setup():
    """Extract tRNAScan covariance models as separate files."""
    base = os.path.join("usr", "local", "lib", "tRNAscan-SE", "models")
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
            parts = line.split("\t")
            data[parts[0].strip()] = {
                "score": float(parts[8].strip()),
                "isotype": parts[4].strip(),
                "anticodon": parts[5].strip(),
                "note": parts[9].lower().strip(),
                "start": int(parts[2].strip()),
                "end": int(parts[3].strip()),
            }
    return data


def run_trnascan(fasta_input, output_folder, domain):
    """Launch tRNAScan-SE and return parsed results."""
    output_file = os.path.join(
        output_folder,
        domain + "-" + os.path.basename(fasta_input).replace(".fasta", ".txt"),
    )
    if domain == "M":
        domain = "M vert"
    if not os.path.exists(output_file):
        cmd = f"tRNAscan-SE -{domain} -o {output_file} {fasta_input}"
        print(cmd)
        os.system(cmd)
    return parse_trnascan_output(output_file)


def skip_trna(entry):
    """
    Some tRNAs should not be drawn and need to be skipped.
    """
    if "pseudo" in entry["note"] or entry["isotype"] in ["Undet", "Sup"]:
        return True
    return False


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


def visualise(
    domain, isotype, fasta_input, output_folder, constraint, exclusion, fold_type
):
    """A wrapper for visualising multiple tRNA sequences in a FASTA file."""
    filename = "headers.txt"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if not os.path.exists(f"{fasta_input}.ssi"):
        cmd = f"esl-sfetch --index {fasta_input}"
        os.system(cmd)

    cmd = f"grep '>' {fasta_input} > {filename}"
    os.system(cmd)

    with open(filename, "r", encoding="utf-8") as f_headers:
        for _, line in enumerate(f_headers):
            seq_id = line.split(" ", 1)[0].replace(">", "").strip()
            print(seq_id)
            ribovision.visualise(
                "gtrnadb",
                fasta_input,
                output_folder,
                seq_id,
                None,
                constraint,
                exclusion,
                fold_type,
                domain,
                isotype,
                None,
                None,
            )
    os.system(f"rm {filename}")


def get_trnascan_cm(domain, isotype):
    """
    Fetch a domain-specific isotype covariance model as a separate file.
    """
    if not os.path.exists(config.GTRNADB_CM_LIBRARY):
        os.mkdir(config.GTRNADB_CM_LIBRARY)
    cm_output = Path(config.GTRNADB_CM_LIBRARY) / f"{domain}_{isotype}.cm"
    if cm_output.exists():
        return str(cm_output)

    cm_library = Path("/usr/local/lib/tRNAscan-SE/models")
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

    with cm_output.open("w") as out:
        cmd = ["cmfetch", str(cm_library), cm_name]
        sp.check_call(cmd, stdout=out)
    return cm_output


def get_traveler_template_xml(domain, isotype):
    """Get Traveler template with coordinates."""
    if domain == "A":
        return os.path.join(
            config.GTRNADB_ARCH, f"arch-{isotype}-traveler-template.xml"
        )
    elif domain == "B":
        return os.path.join(
            config.GTRNADB_BACT, f"bact-{isotype}-traveler-template.xml"
        )
    elif domain == "M":
        if "Leu" in isotype or "Ser" in isotype:
            isotype = isotype[0:3] + "_" + isotype[3:6]
        return os.path.join(
            config.GTRNADB_MITO, f"mito_vert_{isotype}-traveler-template.xml"
        )
    elif domain == "E":
        return os.path.join(config.GTRNADB_EUK, f"euk-{isotype}-traveler-template.xml")
    else:
        raise ValueError(f"Unknown domain {domain}")


def get_traveler_fasta(domain, isotype):
    """Get Traveler structure file."""
    if domain == "A":
        return os.path.join(config.GTRNADB_ARCH, f"arch-{isotype}-traveler.fasta")
    elif domain == "B":
        return os.path.join(config.GTRNADB_BACT, f"bact-{isotype}-traveler.fasta")
    elif domain == "M":
        if "Leu" in isotype or "Ser" in isotype:
            isotype = isotype[0:3] + "_" + isotype[3:6]
        return os.path.join(config.GTRNADB_MITO, f"mito_vert_{isotype}-traveler.fasta")
    elif domain == "E":
        return os.path.join(config.GTRNADB_EUK, f"euk-{isotype}-traveler.fasta")
    else:
        raise ValueError(f"Unknown domain {domain}")
