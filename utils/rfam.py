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

import glob
import gzip
import io
import os
import re
import shutil
import tempfile
from pathlib import Path

import requests
from tqdm import tqdm

from . import config, core
from .rfamseed import RfamSeed
from .ribovore import MIN_GA
from .rnartist import RnaArtist
from .rnartist_setup import compare_rnartist_and_rscape
from .runner import runner
from .scale_template import scale_coordinates
from .shared import cmfetch

# these RNAs are better handled by other methods
BLACKLIST = [
    "RF00001",  # 5S
    "RF00002",  # 5.8S
    "RF02541",  # LSU_rRNA_bacteria
    "RF00177",  # SSU_rRNA_bacteria
    "RF01960",  # SSU_rRNA_eukarya
    "RF02540",  # LSU_rRNA_archaea
    "RF02543",  # LSU_rRNA_eukarya
    "RF01959",  # SSU_rRNA_archaea
    "RF02542",  # SSU_rRNA_microsporidia
    "RF02546",  # LSU_trypano_mito
    "RF02545",  # SSU_trypano_mito
    "RF00005",  # tRNA
    "RF01852",  # tRNA-Sec
    "RF00009",  # Nuclear RNase P
    "RF00010",  # Bacterial RNase P class A
    "RF00011",  # Bacterial RNase P class B
    "RF00373",  # Archaeal RNase P
    "RF01577",  # Plasmodium RNase_P
    "RF00061",  # HCV IRES
    "RF02357",  # RNAse P Type T
    "RF00023",  # tmRNA
    "RF01849",  # Alphaproteobacteria transfer-messenger RNA
    "RF01850",  # Betaproteobacteria transfer-messenger RNA
    "RF01851",  # Cyanobacteria transfer-messenger RNA
    "RF02544",  # Mitochondrion-encoded tmRNA
]

PREFER_RNARTIST_LIST = Path(config.RFAM_DATA) / "prefer_rnartist.txt"


def blacklisted():
    """Get a list of blacklisted Rfam families."""
    bad = set(BLACKLIST)
    filename = Path(config.RFAM_DATA) / "no_structure.txt"
    with open(filename, "r", encoding="utf-8") as raw:
        bad.update(l.strip() for l in raw)
    return bad


def get_traveler_template_xml(rfam_acc, method="auto"):
    """Get a path to a template file given an Rfam accession."""

    if method not in ["auto", "r2r", "rnartist", "rscape"]:
        raise ValueError(f"Unknown method: {method}")

    rfam_data_path = Path(config.RFAM_DATA) / rfam_acc
    traveler_path = rfam_data_path / "traveler-template.xml"
    rnartist_path = rfam_data_path / "rnartist-template.xml"

    if method in ["r2r", "rscape"]:
        return traveler_path
    if method == "rnartist":
        return rnartist_path

    if PREFER_RNARTIST_LIST.exists():
        with PREFER_RNARTIST_LIST.open("r") as f_in:
            if rfam_acc in (line.strip() for line in f_in):
                return rnartist_path
    return traveler_path


def get_traveler_fasta(rfam_acc):
    """Get a path to a consensus structure FASTA given an Rfam accession."""
    filename = os.path.join(config.RFAM_DATA, rfam_acc, f"{rfam_acc}-traveler.fasta")
    return filename


def get_rfam_cm(rfam_acc):
    """Get a path to an Rfam covariance model given an accession."""
    if rfam_acc == "RF00005":
        return str(Path(config.RFAM_DATA) / rfam_acc / f"{rfam_acc}.cm")
    return cmfetch(rfam_acc)


# pylint: disable=too-many-locals, too-many-statements, too-many-branches
def get_rfam_cms():
    """Fetch Rfam covariance models excluding blacklisted models."""
    rfam_all_cm = Path(config.RFAM_CM_LIBRARY) / "all.cm"

    print("Deleting old Rfam library")
    rfam_all_cm.unlink(missing_ok=True)

    print("Downloading Rfam.cm from Rfam FTP")
    rfam_cm = Path(config.RFAM_DATA) / "Rfam.cm"
    rfam_ids = Path(config.RFAM_DATA) / "rfam_ids.txt"
    if not rfam_cm.exists():
        url = "http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
        rfam_gz = rfam_cm.with_suffix(".gz")

        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for HTTP errors

        with rfam_gz.open("wb") as out_file:
            for chunk in response.iter_content(chunk_size=8192):
                out_file.write(chunk)

        with gzip.open(rfam_gz, "rb") as f_in:
            with rfam_gz.with_suffix(".cm").open("wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        rfam_gz.unlink()  # Remove the .gz file after decompression

    rfam_ssi = rfam_cm.with_suffix(rfam_cm.suffix + ".ssi")
    if not rfam_ssi.exists():
        print("Indexing Rfam.cm")
        runner.run(f"cmfetch --index {rfam_cm}")

    if not rfam_ids.exists():
        print("Extracting Rfam IDs")
        with open(rfam_cm, "r") as infile, open(rfam_ids, "w") as outfile:
            for line in infile:
                if line.startswith("ACC   RF"):
                    parts = line.split()
                    if len(parts) > 1:
                        outfile.write(parts[1] + "\n")

    print("Fetching whitelisted Rfam CMs")
    rfam_accs = []
    with open(rfam_ids, "r", encoding="utf-8") as f_in:
        for line in f_in:
            rfam_acc = line.strip()
            if rfam_acc in blacklisted():
                continue
            rfam_accs.append(rfam_acc)

    # Fetch the covariance models
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_list:
        for rfam_acc in rfam_accs:
            temp_list.write(f"{rfam_acc}\n")
        temp_list.flush()
        runner.run(f"cmfetch -f {rfam_cm} {temp_list.name} >> {rfam_all_cm}")

    # Get the accession to Rfam ID mapping
    family_file_path = Path(config.RFAM_DATA) / "family.txt"

    # Download the file if it doesn't exist
    if not family_file_path.exists():
        url = "http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz"
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for HTTP errors

        gzipped_path = family_file_path.with_suffix(".txt.gz")
        with gzipped_path.open("wb") as out_file:
            for chunk in response.iter_content(chunk_size=8192):
                out_file.write(chunk)

        # Decompress the file
        with gzip.open(gzipped_path, "rt", errors="ignore") as f_in:
            with family_file_path.open("w") as f_out:
                for line in f_in:
                    f_out.write(line)

        gzipped_path.unlink()  # Remove the .gz file after decompression

    model_ids = {}
    with family_file_path.open(errors="ignore") as f_db_dump:
        for line in f_db_dump:
            if line.startswith("RF"):
                parts = line.split()
                model_id = parts[1]
                rfam_acc = parts[0]
                model_ids[rfam_acc] = model_id

    # Generate model info file
    modelinfo = Path(config.RFAM_CM_LIBRARY) / "modelinfo.txt"
    with open(modelinfo, "w", encoding="utf-8") as model_out:
        model_out.write("*all*    -    -    all.cm\n")
        for rfam_acc in rfam_accs:
            model_out.write(
                "    ".join([model_ids[rfam_acc], "SSU", "Bacteria", "-\n"])
            )

    rfam_cm.unlink(missing_ok=True)
    rfam_ssi.unlink(missing_ok=True)


def setup_trna_cm():
    """Make sure the RF00005 tRNA model exists as it is used as a fallback
    for all tRNA models that do not match tRNAScan-SE."""
    rfam_acc = "RF00005"
    trna_cm_path = Path(config.RFAM_DATA) / rfam_acc / f"{rfam_acc}.cm"

    # Create the directory if it doesn't exist
    trna_cm_path.parent.mkdir(parents=True, exist_ok=True)

    # Download the file if it doesn't exist
    if not trna_cm_path.exists():
        url = f"https://rfam.org/family/{rfam_acc}/cm"
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for HTTP errors
        trna_cm_path.write_bytes(response.content)

        rscape2traveler(rfam_acc)

        if not trna_cm_path.exists():
            raise FileNotFoundError(f"Rfam tRNA CM not found in {trna_cm_path}")


def delete_preexisting_rfam_data():
    """Delete preexisting Rfam data."""
    # delete Rfam cms
    os.system(f"rm -f {config.RFAM_CM_LIBRARY}/all.cm")
    os.system(f"rm -f {config.RFAM_CM_LIBRARY}/all.cm.ssi")
    os.system(f"rm -f {config.RFAM_CM_LIBRARY}/all.cm.i1*")
    os.system(f"rm -f {config.RFAM_CM_LIBRARY}/modelinfo.txt")
    # delete full Rfam.cm downloaded from FTP
    os.system(f"rm -f {config.RFAM_DATA}/Rfam.cm")
    os.system(f"rm -f {config.RFAM_DATA}/Rfam.cm.ssi")
    # delete Rfam seed
    os.system(f"rm -f {config.RFAM_CM_LIBRARY}/Rfam.seed")
    os.system(f"rm -f {config.RFAM_CM_LIBRARY}/Rfam.seed.ssi")
    # delete tmp cache
    os.system("rm -Rf /tmp/rfam_seed")
    os.system("rm -Rf /tmp/cms")
    # delete template files
    os.system(f"rm -Rf {config.RFAM_DATA}/RF0*")
    # delete summary files
    os.system(f"rm -Rf {config.RFAM_DATA}/family.txt")
    os.system(f"rm -Rf {config.RFAM_DATA}/rfam_ids.txt")


def setup_rnartist(rerun=False):
    """Generate a list of Rfam accessions where RNArtist templates are preferred."""
    # Clear cached seeds and CMs to force re-extraction from the current archives
    os.system("rm -Rf /tmp/rfam_seed")
    os.system("rm -Rf /tmp/cms")
    prefer_rnartist = []
    for rfam_acc in tqdm(get_all_rfam_acc(), "Comparing R-scape and RNArtist"):
        if rfam_acc in BLACKLIST:
            continue
        if not has_structure(rfam_acc):
            continue
        rscape2traveler(rfam_acc)
        RnaArtist(rfam_acc).run(rerun=rerun)
        chosen_template, overlaps = compare_rnartist_and_rscape(rfam_acc)
        if chosen_template == "rnartist":
            prefer_rnartist.append(rfam_acc)
            tqdm.write(
                f"Prefer RNArtist for {rfam_acc}. "
                f"R-scape overlaps: {overlaps['rscape']}, "
                f"RNArtist overlaps: {overlaps['rnartist']}"
            )
    with open(PREFER_RNARTIST_LIST, "w") as f_out:
        for accession in prefer_rnartist:
            f_out.write(f"{accession}\n")


def setup(accessions=None) -> None:
    """Setup Rfam template library."""
    delete_preexisting_rfam_data()
    get_rfam_cms()
    RfamSeed().download_rfam_seed_archive().get_no_structure_file()
    if not accessions:
        accessions = get_all_rfam_acc()
    for accession in tqdm(accessions, "Generating R2R templates"):
        if accession in BLACKLIST:
            continue
        if not has_structure(accession):
            continue
        rscape2traveler(accession)
    setup_trna_cm()
    # delete temporary files
    os.system(f"cd {config.RFAM_DATA} && ./clean_up_files.sh")


# pylint: disable-next=too-many-branches
def generate_traveler_fasta(rfam_acc):
    """
    Generate fasta format for Rfam consensus.

    Example (truncated):

    >RF00162
    NNCUUAUCNAGAGNGGYRGAGGGAYNGGCCC
    ((((((((......(((...(((.....)))
    """
    ss_cons = ""
    consensus = ""

    # get a list of alignments
    seeds = []
    for seed in glob.glob(os.path.join(config.RFAM_DATA, rfam_acc, "*.R2R.sto")):
        seeds.append(seed)
    if len(seeds) != 1:
        print("Error: unusual number of seed alignments")

    with open(seeds[0], "r", encoding="utf-8") as f_seed:
        for line in f_seed.readlines():
            if line.startswith("#=GC SS_cons "):
                parts = line.split()
                ss_cons += parts[2].replace("<", "(").replace(">", ")")
            elif line.startswith("#=GC cons") and "conss" not in line:
                parts = line.split()
                consensus += parts[2]

    if not ss_cons:
        print("No SS_CONS found")
    elif len(ss_cons) == len(consensus):
        if "-" in consensus or "." in consensus:
            new_ss_cons = []
            new_consensus = []
            for i, nucleotide in enumerate(consensus):
                if nucleotide == "-" and ss_cons[i] == ".":
                    pass
                elif nucleotide == "-" and ss_cons[i] in "()":
                    # RF00016 for example
                    new_ss_cons.append(ss_cons[i])
                    new_consensus.append("n")
                else:
                    new_ss_cons.append(ss_cons[i])
                    new_consensus.append(consensus[i])
            ss_cons = "".join(new_ss_cons)
            consensus = "".join(new_consensus)

        with open(get_traveler_fasta(rfam_acc), "w", encoding="utf-8") as f_seed:
            f_seed.write(f">{rfam_acc}\n")
            f_seed.write(f"{consensus.upper()}\n")
            f_seed.write(f"{ss_cons}\n")
    else:
        print("Error: structure and consensus have different lengths")


def convert_path_to_text(line):
    # pylint: disable=line-too-long
    """
    <!--
    <path
     fill="#d90000" stroke="#000000" stroke-width="0.72" stroke-linecap="butt" stroke-linejoin="miter"
     d="M 91.0717 207.452 A 2.5575,2.5575 0 0,1 85.9567,207.452A 2.5575,2.5575 0 0,1 91.0717,207.452Z"/>
     -->

    <text x="85.9567" y="210.352" id="text1002">
      <tspan fill="#807b88"  font-variant="normal" font-weight="normal" font-style="normal" font-family="Bitstream Vera Sans" font-size="7.5" id="tspan1003">N</tspan>
    </text>
    """
    match = re.search(r'd="M (\d+(\.\d+)?) (\d+(\.\d+)?) ', line)
    if match:
        x_coord = float(match.group(1))
        y_coord = float(match.group(3))
        new_x = x_coord - 5.115
        new_y = y_coord + 2.9

        text = """
        <text x="{}" y="{}" id="foobar">
          <tspan fill="#807b88"  font-variant="normal" font-weight="normal" font-style="normal" font-family="Bitstream Vera Sans" font-size="7.5">N</tspan>
        </text>"""

        xml = """<point x="{:.2f}" y="{:.2f}" b="N"/>\n"""

        return (text.format(new_x, new_y), xml.format(new_x, new_y))
    print(line)
    print("convert_path_to_text did not find a match")
    return ""


def convert_text_to_xml(line):
    # pylint: disable=line-too-long
    """
    # <text x="209.519" y="231.111" id="text1002"><tspan x="209.519" y="231.111" fill="#807b88"  font-variant="normal" font-weight="normal" font-style="normal" font-family="Bitstream Vera Sans" font-size="7.5" id="tspan1003">R</tspan></text>
    # <point x="209.52" y="0.52" b="231.111"/>

    <text x="85.8067" y="195.025" id="text1002"><tspan x="85.8067" y="195.025" fill="#807b88"  font-variant="normal" font-weight="normal" font-style="normal" font-family="Bitstream Vera Sans" font-size="7.5" id="tspan1003">C</tspan></text>
    """
    match = re.search(
        r'<text x="(\d+(\.\d+)?)" y="(\d+(\.\d+)?)".+>(.+?)</tspan></text>', line
    )
    if match:
        point = '<point x="{:.2f}" y="{:.2f}" b="{}"/>\n'
        return point.format(
            float(match.group(1)), float(match.group(3)), match.group(5)
        )
    print(line)
    print("convert_text_to_xml did not find a match")
    return ""


def get_all_rfam_acc():
    """Get a list of Rfam accessions from an FTP database dump file."""
    rfam_accs = []
    family_file_path = Path(config.RFAM_DATA) / "family.txt"

    # Download the file if it doesn't exist
    if not family_file_path.exists():
        url = "http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz"
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for HTTP errors

        gzipped_path = family_file_path.with_suffix(".txt.gz")
        with gzipped_path.open("wb") as out_file:
            for chunk in response.iter_content(chunk_size=8192):
                out_file.write(chunk)

        # Decompress the file
        with gzip.open(gzipped_path, "rt", errors="ignore") as f_in:
            with family_file_path.open("w") as f_out:
                for line in f_in:
                    f_out.write(line)

        gzipped_path.unlink()  # Remove the .gz file after decompression

    with family_file_path.open(errors="ignore") as f_db_dump:
        for line in f_db_dump:
            if line.startswith("RF"):
                rfam_acc = line[:7]
                if (
                    rfam_acc in BLACKLIST
                ):  # Assuming BLACKLIST is defined somewhere in your code
                    continue
                rfam_accs.append(rfam_acc)
    return rfam_accs


def get_rfam_acc_by_id(rfam_id):
    """Get Rfam accession corresponding to an Rfam ID.
    Example: return RF00162 for SAM riboswitch."""
    family_file_path = Path(config.RFAM_DATA) / "family.txt"
    # Download and decompress the file if it doesn't exist
    if not family_file_path.exists():
        url = "http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz"
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for HTTP errors

        gzipped_path = family_file_path.with_suffix(".txt.gz")
        with gzipped_path.open("wb") as out_file:
            for chunk in response.iter_content(chunk_size=8192):
                out_file.write(chunk)

        # Decompress the file
        with gzip.open(gzipped_path, "rt", errors="ignore") as f_in:
            with family_file_path.open("w") as f_out:
                for line in f_in:
                    f_out.write(line)

        gzipped_path.unlink()  # Remove the .gz file after decompression

    with family_file_path.open(errors="ignore") as raw:
        for line in raw:
            parts = line.split()
            if parts[1] == rfam_id:
                return parts[0]

    raise ValueError(f"Cannot find Rfam accession for: {rfam_id}")


def remove_pseudoknot_from_ss_cons(rfam_seed):
    """
    Some Rfam alignments contain pseudoknot annotations in SS_cons. To avoid
    R-scape displaying them on diagrams, pseudoknots need to be removed before
    running R-scape.
    """
    seed_no_pk = rfam_seed.with_suffix(".nopk.sto")
    with io.open(rfam_seed, "r", encoding="latin-1") as f_seed_in:
        with open(seed_no_pk, "w", encoding="utf-8") as f_seed_out:
            for line in f_seed_in.readlines():
                if line.startswith("#=GC SS_cons "):
                    match = re.match(r"(#=GC SS_cons)(\s+)(.+)", line)
                    no_pk = re.sub(r"\w", ".", match.group(3))
                    f_seed_out.write(
                        "".join([match.group(1), match.group(2), no_pk, "\n"])
                    )
                else:
                    f_seed_out.write(line)
    return seed_no_pk


def run_rscape(rfam_acc, destination):
    """
    Run R-scape on Rfam seed alignment to get the R-scape/R2R layout.
    """
    rfam_seed = RfamSeed().get_rfam_seed(rfam_acc)
    rfam_seed_no_pk = remove_pseudoknot_from_ss_cons(rfam_seed)
    if not os.path.exists(os.path.join(destination, "rscape.done")):
        cmd = "R-scape --outdir {folder} {rfam_seed} && touch {folder}/rscape.done".format(
            folder=destination, rfam_seed=rfam_seed_no_pk
        )
        runner.run(cmd)
        # delete any temporary r2r_meta files
        for filename in glob.glob("*.r2r_meta"):
            os.remove(filename)

    rscape_svg = None
    for svg in glob.glob(os.path.join(destination, "*.svg")):
        if "R2R.sto.svg" in svg:
            rscape_svg = svg
    if not rscape_svg:
        print("Error: R-scape SVG file not found")
    return rscape_svg


def convert_rscape_svg_to_one_line(rscape_svg, destination):
    """
    Convert R-scape SVG into SVG with 1 line per element.
    """
    # pylint: disable=consider-using-f-string
    output = os.path.join(destination, "rscape-one-line.svg")
    cmd = (
        r"perl -0777 -pe 's/\n +fill/ fill/g' {rscape_svg} | "
        r"perl -0777 -pe 's/\n d=/ d=/g' | "
        r"perl -0777 -pe 's/\n +<tspan/ <tspan/g' | "
        r"perl -0777 -pe 's/\n<\/text>/<\/text>/g' "
        r"> {output}"
    ).format(rscape_svg=rscape_svg, output=output)
    runner.run(cmd)
    return output


# pylint: disable=too-many-branches
def convert_rscape_svg_to_traveler(rscape_one_line_svg, destination):
    """
    Convert R-scape SVG into traveler xml and SVG.
    """
    header = """
    <svg
    xmlns:svg="http://www.w3.org/2000/svg"
    xmlns="http://www.w3.org/2000/svg"
    xmlns:xlink="http://www.w3.org/1999/xlink"
    version="1.0"
    width="1000"
    height="2000"
    xml:space="preserve">
    """
    footer = "</svg>"

    xml_header = "<structure>\n"
    xml_footer = "</structure>"

    traveler_template_svg = os.path.join(destination, "traveler-template.svg")
    with open(rscape_one_line_svg, "r", encoding="utf-8") as f_in:
        with open(traveler_template_svg, "w", encoding="utf-8") as f_out:
            with open(
                os.path.join(destination, "traveler-template.xml"),
                "w",
                encoding="utf-8",
            ) as xml_out:
                f_out.write(header)
                xml_out.write(xml_header)

                for line in f_in.readlines():
                    if "#5c5c5c" in line:
                        continue
                    if "text1000" in line:
                        continue
                    if "#d7efc5" in line:
                        continue
                    if "#31a354" in line:  # lancaster helix covariation
                        continue
                    if '<path fill="none" stroke="#000000"' in line:
                        continue
                    # circles indicating GU pairs
                    if '<path fill="#000000" stroke="none"' in line:
                        continue
                    if "&apos;" in line:
                        continue
                    if "5&#x2032" in line:  # 5' label
                        continue
                    if "pk" in line:
                        continue
                    if line.startswith("<path"):
                        text, xml = convert_path_to_text(line)
                        f_out.write(text)
                        xml_out.write(xml)
                    elif line.startswith("<text"):
                        xml = convert_text_to_xml(line)
                        if not xml:
                            continue
                        xml_out.write(xml)
                        f_out.write(line)
                    else:
                        continue
                f_out.write(footer)
                xml_out.write(xml_footer)


def rscape2traveler(rfam_acc):
    """Create a Traveler XML template based on an
    R-scape consensus 2D layout."""
    destination = os.path.join(config.RFAM_DATA, rfam_acc)
    if not os.path.exists(destination):
        os.makedirs(destination)
    if os.path.exists(get_traveler_fasta(rfam_acc)) and os.path.exists(
        get_traveler_template_xml(rfam_acc, "r2r")
    ):
        return
    rscape_svg = run_rscape(rfam_acc, destination)
    rscape_one_line_svg = convert_rscape_svg_to_one_line(rscape_svg, destination)
    convert_rscape_svg_to_traveler(rscape_one_line_svg, destination)
    scale_coordinates(get_traveler_template_xml(rfam_acc, "r2r"))
    generate_traveler_fasta(rfam_acc)


# pylint: disable-next=too-many-arguments,too-many-positional-arguments
def generate_2d(
    rfam_acc,
    output_folder,
    fasta_input,
    constraint,
    exclusion,
    fold_type,
    quiet=False,
    rfam_template_type="rscape",
):
    """Loop over the sequences in fasta file and visualise each
    using the family template."""
    destination = f"{output_folder}/{rfam_acc}"
    if not os.path.exists(destination):
        os.makedirs(destination)

    if not os.path.exists(fasta_input + ".ssi"):
        runner.run(f"esl-sfetch --index {fasta_input}")
    # pylint: disable=consider-using-with
    headers = tempfile.NamedTemporaryFile(delete=False).name
    with open(fasta_input, "r") as infile, open(headers, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                outfile.write(line)

    with open(headers) as f_headers:
        for line in f_headers:
            seq_id = line.split(" ", 1)[0].replace(">", "").strip()
            core.visualise(
                "rfam",
                fasta_input,
                destination,
                seq_id,
                rfam_acc,
                constraint,
                exclusion,
                fold_type,
                domain=None,
                isotype=None,
                start=None,
                end=None,
                quiet=quiet,
                rfam_template=rfam_template_type,
            )
    Path(headers).unlink(missing_ok=True)


def has_structure(rfam_acc):
    """Return a list of families that have consensus 2D structure."""
    no_structure = []
    no_structure_filename = os.path.join(config.RFAM_DATA, "no_structure.txt")
    with open(no_structure_filename) as f_list:
        for line in f_list.readlines():
            no_structure.append(line.strip())
    return rfam_acc not in no_structure


def cmsearch_nohmm_mode(fasta_input, output_folder, rfam_acc):
    """
    Run cmsearch on the fasta sequences using cmsearch in the --nohmm mode
    to get potentially missing hits.
    """
    subfolder = os.path.join(output_folder, rfam_acc)
    os.makedirs(subfolder, exist_ok=True)
    tblout = os.path.join(subfolder, "cmsearch.tblout")
    outfile = os.path.join(subfolder, "cmsearch.out.txt")
    cm_file = os.path.join(config.RFAM_DATA, rfam_acc, f"{rfam_acc}.cm")
    runner.run(
        f"cmsearch --nohmm -o {outfile} --tblout {tblout} {cm_file} {fasta_input}"
    )
    hits = os.path.join(subfolder, "hits.txt")
    with open(tblout, "r") as infile, open(hits, "w") as outfile:
        for line in infile:
            if not line.startswith("#") and "?" not in line:
                parts = line.split()
                if len(parts) >= 14:
                    bit_score = float(parts[14])
                    strand = parts[9]
                    if bit_score >= MIN_GA and strand == "+":
                        outfile.write(f"{parts[0]}\t{parts[3]}\tPASS\n")
    ids = set()
    with open(hits) as f_hits:
        for line in f_hits:
            if "\t" in line:
                hit_id, _, _ = line.strip().split("\t")
                ids.add(hit_id)
    return ids
