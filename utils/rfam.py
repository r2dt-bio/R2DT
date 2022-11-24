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
import io
import os
import re
import subprocess as sp

from . import config, core
from . import generate_model_info as mi

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
    "RF00061",  # HCV IRES
    "RF02357",  # RNAse P Type T
]


def blacklisted():
    """Get a list of blacklisted Rfam families."""
    bad = set(BLACKLIST)
    filename = os.path.join(config.RFAM_DATA, "no_structure.txt")
    with open(filename, "r", encoding="utf-8") as raw:
        bad.update(l.strip() for l in raw)
    return bad


def get_traveler_template_xml(rfam_acc):
    """Get a path to a template file given an Rfam accession."""
    filename = os.path.join(config.RFAM_DATA, rfam_acc, "traveler-template.xml")
    return filename


def get_traveler_fasta(rfam_acc):
    """Get a path to a consensus structure FASTA given an Rfam accession."""
    filename = os.path.join(config.RFAM_DATA, rfam_acc, f"{rfam_acc}-traveler.fasta")
    return filename


def get_rfam_cm(rfam_acc):
    """Get a path to an Rfam covariance model given an accession."""
    if rfam_acc == "RF00005":
        return os.path.join(config.RFAM_DATA, rfam_acc, rfam_acc + ".cm")
    return os.path.join(config.CM_LIBRARY, "rfam", rfam_acc + ".cm")


def get_rfam_cms():
    """Fetch Rfam covariance models excluding blacklisted models."""
    rfam_cm_location = os.path.join(config.CM_LIBRARY, "rfam")
    rfam_whitelisted_cm = os.path.join(rfam_cm_location, "all.cm")
    print("Deleting old Rfam library")
    cmd = f"rm -Rf {rfam_cm_location} && mkdir {rfam_cm_location}"
    os.system(cmd)
    print("Downloading Rfam.cm from Rfam FTP")
    rfam_cm = os.path.join(config.RFAM_DATA, "Rfam.cm")
    rfam_ids = os.path.join(config.RFAM_DATA, "rfam_ids.txt")
    if not os.path.exists(rfam_cm):
        cmd = (
            f"wget -O {rfam_cm}.gz "
            f"ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz "
            f"&& gunzip {rfam_cm}.gz"
        )
        os.system(cmd)
    print("Indexing Rfam.cm")
    if not os.path.exists(f"{rfam_cm}.ssi"):
        os.system(f"cmfetch --index {rfam_cm}")
    print("Get a list of all Rfam ids")
    if not os.path.exists(rfam_ids):
        cmd = f"awk '/ACC   RF/ {{print $2}}' {rfam_cm} > {rfam_ids}"
        os.system(cmd)
    print("Fetching whitelisted Rfam CMs")
    with open(rfam_ids, "r", encoding="utf-8") as f_in:
        for line in f_in:
            rfam_acc = line.strip()
            if rfam_acc in blacklisted():
                continue
            print(rfam_acc)
            cmd = f"cmfetch {rfam_cm} {rfam_acc} >> {rfam_whitelisted_cm}"
            os.system(cmd)
            cm_file = os.path.join(rfam_cm_location, f"{rfam_acc}.cm")
            cmd = f"cmfetch {rfam_cm} {rfam_acc} > {cm_file}"
            os.system(cmd)
    print("Cleaning up")
    os.system(f"rm {rfam_cm}")
    os.system(f"rm {rfam_cm}.ssi")


def setup_trna_cm():
    """Make sure the RF00005 tRNA model exists as it is used as a fallback
    for all tRNA models that do not match tRNAScan-SE."""
    rfam_acc = "RF00005"
    trna_cm = os.path.join(config.RFAM_DATA, rfam_acc, f"{rfam_acc}.cm")
    os.system(f"mkdir -p {os.path.join(config.RFAM_DATA, rfam_acc)}")
    if not os.path.exists(trna_cm):
        cmd = f"wget -O {trna_cm} https://rfam.org/family/{rfam_acc}/cm"
        os.system(cmd)
        rscape2traveler(rfam_acc)
        if not os.path.exists(trna_cm):
            raise Exception(f"Rfam tRNA CM not found in {trna_cm}")


def setup(accessions=None):
    """Setup Rfam template library."""
    get_rfam_cms()
    mi.generate_model_info(cm_library=os.path.join(config.CM_LIBRARY, "rfam"))
    if not accessions:
        accessions = get_all_rfam_acc()
    for accession in accessions:
        if accession in BLACKLIST:
            continue
        rscape2traveler(accession)
    setup_trna_cm()


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
        r'<text x="(\d+(\.\d+)?)" y="(\d+(\.\d+)?)".+>(\w)</tspan></text>', line
    )
    if match:
        point = '<point x="{:.2f}" y="{:.2f}" b="{}"/>\n'
        return point.format(
            float(match.group(1)), float(match.group(3)), match.group(5)
        )
    print(line)
    print("convert_text_to_xml did not find a match")
    return ""


def download_rfam_seed(rfam_acc):
    """Fetch Rfam seed alignment using the API."""
    output = os.path.join(config.RFAM_DATA, rfam_acc, f"{rfam_acc}.seed")
    if not os.path.exists(output):
        url = f"https://rfam.org/family/{rfam_acc}/alignment"
        cmd = f"wget -O {output} {url}"
        os.system(cmd)
    return output


def get_all_rfam_acc():
    """Get a list of Rfam accessions from an FTP database dump file."""
    rfam_accs = []
    family_file = os.path.join(config.RFAM_DATA, "family.txt")
    if not os.path.exists(family_file):
        cmd = (
            f"wget -O {family_file}.gz "
            f"ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz"
        )
        try:
            sp.check_output([cmd], shell=True, stderr=sp.STDOUT)
        except sp.CalledProcessError as error:
            print(f"Error {error.output}")

        if not os.path.exists(f"{family_file}.gz"):
            print("Family file not downloaded")

        cmd = f"gunzip {family_file}.gz"
        try:
            sp.check_output([cmd], shell=True, stderr=sp.STDOUT)
        except sp.CalledProcessError as error:
            print(f"Error {error.output}")
    with open(family_file, encoding="utf8", errors="ignore") as f_db_dump:
        for line in f_db_dump:
            if line.startswith("RF"):
                rfam_acc = line[:7]
                if rfam_acc in BLACKLIST:
                    continue
                rfam_accs.append(rfam_acc)
    print(f"Found {len(rfam_accs)} Rfam accessions")
    return rfam_accs


def get_rfam_acc_by_id(rfam_id):
    """Get Rfam accession corresponding to an Rfam ID.
    Example: return RF00162 for SAM riboswitch."""
    family_file = os.path.join(config.RFAM_DATA, "family.txt")
    if not os.path.exists(family_file):
        cmd = (
            f"wget -O {family_file}.gz "
            f"ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz && "
            f"gunzip {family_file}.gz"
        )
        os.system(cmd)

    with open(family_file, encoding="utf8", errors="ignore") as raw:
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
    seed_no_pk = os.path.join(
        os.path.dirname(rfam_seed), f"nopk-{os.path.basename(rfam_seed)}"
    )
    with io.open(rfam_seed, "r", encoding="latin-1") as f_seed_in:
        with open(seed_no_pk, "w", encoding="utf-8") as f_seed_out:
            for line in f_seed_in.readlines():
                if line.startswith("#=GC SS_cons "):
                    match = re.match(r"(#=GC SS_cons)(\s+)(.+)", line)
                    no_pk = re.sub(r"\w", ".", match.group(3))
                    # import pdb; pdb.set_trace()
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
    rfam_seed = download_rfam_seed(rfam_acc)
    rfam_seed_no_pk = remove_pseudoknot_from_ss_cons(rfam_seed)
    if not os.path.exists(os.path.join(destination, "rscape.done")):
        cmd = "R-scape --outdir {folder} {rfam_seed} && touch {folder}/rscape.done".format(
            folder=destination, rfam_seed=rfam_seed_no_pk
        )
        os.system(cmd)

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
    os.system(cmd)
    return output


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

    with open(rscape_one_line_svg, "r", encoding="utf-8") as f_in:
        with open(
            os.path.join(destination, "traveler-template.svg"), "w", encoding="utf-8"
        ) as f_out:
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
                    if '<path fill="none" stroke="#000000"' in line:
                        continue
                    if "&apos;" in line:
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
        get_traveler_template_xml(rfam_acc)
    ):
        return
    rscape_svg = run_rscape(rfam_acc, destination)
    rscape_one_line_svg = convert_rscape_svg_to_one_line(rscape_svg, destination)
    convert_rscape_svg_to_traveler(rscape_one_line_svg, destination)
    generate_traveler_fasta(rfam_acc)


# pylint: disable-next=too-many-arguments
def generate_2d(rfam_acc, output_folder, fasta_input, constraint, exclusion, fold_type):
    """Loop over the sequences in fasta file and visualise each
    using the family template."""
    destination = f"{output_folder}/{rfam_acc}"
    if not os.path.exists(destination):
        os.makedirs(destination)

    if not os.path.exists(fasta_input + ".ssi"):
        cmd = f"esl-sfetch --index {fasta_input}"
        os.system(cmd)

    headers = "headers.txt"
    cmd = f"grep '>' {fasta_input} > {headers}"
    os.system(cmd)

    with open(headers, "r", encoding="utf-8") as f_headers:
        for line in f_headers:
            seq_id = line.split(" ", 1)[0].replace(">", "").strip()
            print(seq_id)
            core.visualise(
                "rfam",
                fasta_input,
                destination,
                seq_id,
                rfam_acc,
                constraint,
                exclusion,
                fold_type,
            )
    os.system(f"rm {headers}")


def has_structure(rfam_acc):
    """Return a list of families that have consensus 2D structure."""
    no_structure = []
    no_structure_filename = os.path.join(config.RFAM_DATA, "no_structure.txt")
    with open(no_structure_filename, "r", encoding="utf-8") as f_list:
        for line in f_list.readlines():
            no_structure.append(line.strip())
    return rfam_acc not in no_structure


def cmsearch_nohmm_mode(fasta_input, output_folder, rfam_acc):
    """
    Run cmsearch on the fasta sequences using cmsearch in the --nohmm mode
    to get potentially missing hits.
    """
    subfolder = os.path.join(output_folder, rfam_acc)
    os.system(f"mkdir -p {subfolder}")
    tblout = os.path.join(subfolder, "cmsearch.tblout")
    outfile = os.path.join(subfolder, "cmsearch.out.txt")
    cm_file = os.path.join(config.RFAM_DATA, rfam_acc, f"{rfam_acc}.cm")
    cmd = f"cmsearch --nohmm -o {outfile} --tblout {tblout} {cm_file} {fasta_input}"
    print(cmd)
    os.system(cmd)
    hits = os.path.join(subfolder, "hits.txt")
    cmd = (
        f"cat {tblout} | grep -v '^#' | grep -v '?' | "
        f"awk -v OFS='\t' '{{print $1, $4, \"PASS\"}}' > {hits}"
    )
    os.system(cmd)
    ids = set()
    with open(hits, "r", encoding="utf-8") as f_hits:
        for line in f_hits:
            if "\t" in line:
                hit_id, _, _ = line.strip().split("\t")
                ids.add(hit_id)
    return ids
