#!/usr/bin/env python3

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

# pylint: disable=too-many-lines
import glob
import json
import os
import re
import shutil
import subprocess
import tarfile
import time
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

import click  # pylint: disable=import-error
from rich import print as rprint

from tests import tests
from utils import cif_basepairs, config, core
from utils import fr3d as fr3d_utils
from utils import generate_cm_library as gcl
from utils import generate_model_info as gmi
from utils import gtrnadb, lbn_export
from utils import list_models as lm
from utils import pdb_fetch, pdb_post, r2r, rfam
from utils import rna2djsonschema as r2djs
from utils import rnapuzzler
from utils import rnaview as rnaview_utils
from utils import shared
from utils import stockholm as stockholm_utils
from utils import viewer_export, viewer_html
from utils import workstation as workstation_mod
from utils.rnartist import RnaArtist
from utils.runner import runner
from utils.scale_template import scale_coordinates


class Timer:
    """
    Context manager that logs execution time.
    """

    def __init__(self, msg: str, quiet: bool = False):
        self.msg = msg
        self.quiet = quiet
        self.start = None
        self.end = None
        self.interval = None

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end = time.time()
        self.interval = self.end - self.start
        if not self.quiet:
            rprint(
                f"[yellow]Elapsed time for {self.msg}[/yellow]: {self.interval:.2f} seconds"
            )


@click.group()
def cli():
    """Required click stub function."""


@cli.command()
def version():
    """
    Print R2DT version information.
    """
    rprint(shared.get_r2dt_version_header())


@cli.command()
def setup():
    """
    Generate all templates from scratch.
    """
    rprint(shared.get_r2dt_version_header())
    crw_setup()
    rfam.setup()
    gtrnadb.setup()


def crw_setup():
    """Setup CRW CM library."""
    if os.path.exists(config.CRW_CM_LIBRARY):
        rprint("Deleting old CRW library")
        shutil.rmtree(config.CRW_CM_LIBRARY)

    # Extract the tar.gz file
    rprint("Extracting precomputed CRW archive")
    with tarfile.open(os.path.join(config.DATA, "crw-cms.tar.gz"), "r:gz") as tar:
        tar.extractall(path=config.DATA)

    # Move the directory
    source_dir = os.path.join(config.DATA, "crw-cms")

    if os.path.exists(source_dir):
        shutil.move(source_dir, config.CRW_CM_LIBRARY)

    # read CRW blacklist
    crw_blacklist = []
    with open(os.path.join(config.DATA, "crw-blacklist.txt")) as f_in:
        for line in f_in:
            if line.startswith("#"):
                continue
            crw_blacklist.append(line.strip())

    # Delete models from the blacklist
    for model in crw_blacklist:
        model_file = os.path.join(config.CRW_CM_LIBRARY, model + ".cm")
        if os.path.exists(model_file):
            os.remove(model_file)

    rprint("Generating CRW modelinfo file")
    gmi.generate_model_info(cm_library=config.CRW_CM_LIBRARY)


@cli.command()
@click.option("--rnartist", default=False, is_flag=True)
def setup_rfam(rnartist):
    """
    Re-generate Rfam templates from scratch.
    """
    rprint(shared.get_r2dt_version_header())
    if not rnartist:
        rprint("Generating Rfam templates")
        rfam.setup()
    rprint("Setting up RNArtist")
    rfam.setup_rnartist(rerun=False)


def get_seq_ids(input_fasta):
    """
    Get a list of sequence ids from a fasta file.
    """
    seq_ids = set()
    with open(input_fasta) as f_in:
        for line in f_in:
            if line.startswith(">"):
                match = re.search(r">(.*?)\s", line)
                if match:
                    seq_ids.add(match.group(1))
    return seq_ids


def get_hits(folder):
    """
    Get a list of sequence ids found in the hits.txt file by ribovore.
    """
    hits = set()
    hits_file = os.path.join(folder, "hits.txt")
    if not os.path.exists(hits_file):
        return hits
    with open(hits_file) as f_in:
        for line in f_in:
            hits.add(line.split("\t")[0])
    return hits


def get_subset_fasta(fasta_input, output_filename, seq_ids):
    """
    Extract a fasta file named <output_filename> with sequence ids <seq_ids>
    from <fasta_input>.
    """
    index_filename = output_filename + ".txt"
    with open(index_filename, "w") as f_out:
        for seq_id in seq_ids:
            f_out.write(f"{seq_id}\n")
    runner.run(f"esl-sfetch -o {output_filename} -f {fasta_input} {index_filename}")
    if not os.path.exists(f"{output_filename}.ssi"):
        runner.run(f"esl-sfetch --index {output_filename}")
    os.remove(index_filename)


def is_templatefree(fasta_input):
    """Check if the input file is a valid fasta file
    with an additional line specifying secondary structure
    in dot bracket format (pseudoknots allowed)."""
    with open(fasta_input) as f_in:
        lines = [line.strip() for line in f_in.readlines() if line.strip()]
    if len(lines) != 3:
        return False
    header, sequence, structure = lines
    if not header.startswith(">"):
        return False
    if len(sequence) != len(structure):
        return False
    if not re.match(r"^[.()<>{}[\]A-z]+$", structure):
        return False
    return True


@cli.command()
@click.argument("fasta-input", type=click.Path())
@click.argument("output-folder", type=click.Path())
@click.option(
    "--force_template",
    type=click.STRING,
    default=None,
    help="Force sequences into a specific template",
)
@click.option(
    "--constraint", default=False, is_flag=True, help="Fold insertions using RNAfold"
)
@click.option("--exclusion", default=None)
@click.option("--fold_type", default=None)
@click.option("--quiet", "-q", default=False, is_flag=True)
@click.option(
    "--skip_ribovore_filters",
    default=False,
    is_flag=True,
    help="Ignore ribovore QC checks",
)
@click.pass_context
def draw(
    ctx,
    fasta_input,
    output_folder,
    force_template,
    constraint,
    exclusion,
    fold_type,
    quiet,
    skip_ribovore_filters,
):
    """
    Single entry point for visualising 2D for an RNA sequence.
    Selects a template and runs Traveler using CRW, LSU, or Rfam libraries.
    """
    # pylint: disable=too-many-arguments,too-many-positional-arguments
    # pylint: disable=too-many-locals,too-many-statements,too-many-branches
    if not quiet:
        rprint(shared.get_r2dt_version_header())

    fasta_input = shared.sanitise_fasta(fasta_input)

    if is_templatefree(fasta_input):
        if not quiet:
            rprint("Detected templatefree input.")
        ctx.invoke(
            templatefree,
            fasta_input=fasta_input,
            output_folder=output_folder,
            quiet=quiet,
        )
        return

    all_seq_ids = get_seq_ids(fasta_input)

    if force_template:
        for seq_id in all_seq_ids:
            force_draw(
                force_template,
                fasta_input,
                output_folder,
                seq_id,
                constraint,
                exclusion,
                fold_type,
                quiet=True,
            )
        return

    os.makedirs(output_folder, exist_ok=True)
    crw_output = os.path.join(output_folder, "crw")
    ribovision_ssu_output = os.path.join(output_folder, "ribovision-ssu")
    ribovision_lsu_output = os.path.join(output_folder, "ribovision-lsu")
    rfam_output = os.path.join(output_folder, "rfam")
    gtrnadb_output = os.path.join(output_folder, "gtrnadb")
    rfam_trna_output = os.path.join(output_folder, "RF00005")
    rnasep_output = os.path.join(output_folder, "rnasep")
    tmrna_output = os.path.join(output_folder, "tmrna")

    hits = set()
    subset_fasta = os.path.join(output_folder, "subset.fasta")
    if not os.path.exists(f"{fasta_input}.ssi"):
        runner.run(f"esl-sfetch --index {fasta_input}")

    def get_output_subfolder(method_name):
        """Get folder within the output folder for a given method."""
        subfolders = {
            "ribovision_draw_ssu": os.path.join(output_folder, "ribovision-ssu"),
            "ribovision_draw_lsu": os.path.join(output_folder, "ribovision-lsu"),
            "rrna_draw": os.path.join(output_folder, "crw"),
            "rnasep_draw": os.path.join(output_folder, "rnasep"),
            "tmrna_draw": os.path.join(output_folder, "tmrna"),
        }
        return subfolders.get(str(method_name), "")

    method_list = [
        "rnasep_draw",
        "tmrna_draw",
        "ribovision_draw_ssu",
        "ribovision_draw_lsu",
        "rrna_draw",  # CRW
    ]
    prev_output_subfolder = None
    for method_name in method_list:
        if prev_output_subfolder:
            hits = hits.union(get_hits(prev_output_subfolder))
            subset = all_seq_ids.difference(hits)
            if subset:
                get_subset_fasta(fasta_input, subset_fasta, subset)
        else:
            subset = all_seq_ids
            shutil.copy(fasta_input, subset_fasta)
            if not os.path.exists(f"{subset_fasta}.ssi"):
                runner.run(f"esl-sfetch --index {subset_fasta}")
        if subset:
            with Timer(f"{method_name}", quiet):
                if not quiet:
                    rprint(f"Analysing {len(subset)} sequences with {method_name}")
                output_subfolder = get_output_subfolder(method_name)
                ctx.invoke(
                    globals()[method_name],
                    fasta_input=subset_fasta,
                    output_folder=output_subfolder,
                    constraint=constraint,
                    exclusion=exclusion,
                    fold_type=fold_type,
                    quiet=True,
                    skip_ribovore_filters=skip_ribovore_filters,
                )
                prev_output_subfolder = output_subfolder

    # Rfam
    hits = hits.union(get_hits(prev_output_subfolder))
    subset = all_seq_ids.difference(hits)
    if not quiet:
        rprint(f"Analysing {len(subset)} sequences with Rfam")
    if subset:
        with Timer("Rfam", quiet):
            with open(
                shared.get_ribotyper_output(
                    subset_fasta,
                    rfam_output,
                    config.RFAM_CM_LIBRARY,
                    skip_ribovore_filters,
                ),
            ) as f_ribotyper:
                for line in f_ribotyper.readlines():
                    seq_id, model_id, _ = line.split("\t")
                    core.visualise(
                        "rfam",
                        subset_fasta,
                        rfam_output,
                        seq_id,
                        model_id,
                        constraint,
                        exclusion,
                        fold_type,
                        domain=None,
                        isotype=None,
                        start=None,
                        end=None,
                        quiet=quiet,
                    )

    # GtRNAdb
    hits = hits.union(get_hits(rfam_output))
    subset = all_seq_ids.difference(hits)
    if subset:
        get_subset_fasta(fasta_input, subset_fasta, subset)
        with Timer("GtRNAdb", quiet):
            if not quiet:
                rprint(f"Analysing {len(subset)} sequences with GtRNAdb")
            for trna in gtrnadb.classify_trna_sequences(subset_fasta, gtrnadb_output):
                core.visualise(
                    "gtrnadb",
                    fasta_input,
                    output_folder + "/gtrnadb",
                    trna["id"],
                    None,
                    constraint,
                    exclusion,
                    fold_type,
                    trna["domain"],
                    trna["isotype"],
                    trna["start"],
                    trna["end"],
                    quiet,
                )

    # Rfam tRNA
    hits = hits.union(get_hits(gtrnadb_output))
    subset = all_seq_ids.difference(hits)
    if subset:
        get_subset_fasta(fasta_input, subset_fasta, subset)
        with Timer("Rfam tRNA", quiet):
            if not quiet:
                rprint(f"Analysing {len(subset)} sequences with Rfam tRNA")
            trna_ids = rfam.cmsearch_nohmm_mode(subset_fasta, output_folder, "RF00005")
            if trna_ids:
                get_subset_fasta(fasta_input, subset_fasta, trna_ids)
                rfam.generate_2d(
                    "RF00005",
                    output_folder,
                    subset_fasta,
                    constraint,
                    exclusion,
                    fold_type,
                    quiet,
                )

    # move svg files to the final location
    result_folders = [
        crw_output,
        ribovision_ssu_output,
        ribovision_lsu_output,
        rfam_output,
        gtrnadb_output,
        rfam_trna_output,
        rnasep_output,
        tmrna_output,
    ]
    for folder in result_folders:
        organise_results(folder, output_folder)
    organise_metadata(output_folder, result_folders)

    # clean up
    os.system(f"rm {output_folder}/subset*")
    os.system(f"rm -f {fasta_input}.ssi")


@cli.command()
def compress_rfam_crw():
    """Generate compressed tar.gz files for the CRW and Rfam all.cm files.
    the files are located in the config.CRW_CM_LIBRARY and config.RFAM_CM_LIBRARY.
    Upon uncompressing the tar.gz files, the files should be also called all.cm.
    I want to only compress the all.cm files, not the entire folder.
    """
    rprint(shared.get_r2dt_version_header())
    rprint("Compressing CRW and Rfam all.cm files")
    crw_cm = Path(config.CRW_CM_LIBRARY) / "all.cm"
    rfam_cm = Path(config.RFAM_CM_LIBRARY) / "all.cm"
    crw_tar = Path(config.CRW_CM_LIBRARY) / "all.cm.tar.gz"
    rfam_tar = Path(config.RFAM_CM_LIBRARY) / "all.cm.tar.gz"
    runner.run(
        f"tar -czf {crw_tar} -C {os.path.dirname(crw_cm)} {os.path.basename(crw_cm)}"
    )
    runner.run(
        f"tar -czf {rfam_tar} -C {os.path.dirname(rfam_cm)} {os.path.basename(rfam_cm)}"
    )
    rprint("Done")


def organise_results(results_folder, output_folder):
    """Move files to the final folder structure."""
    folders = {}
    labels = ["svg", "fasta", "json", "thumbnail"]
    destination = os.path.join(output_folder, "results")
    for label in labels:
        folders[label] = os.path.join(destination, label)
        os.makedirs(folders[label], exist_ok=True)
    svgs = glob.glob(os.path.join(results_folder, "*.svg"))
    if not svgs:
        return
    for svg in svgs:
        if "colored" not in svg:
            continue
        if "enriched" in svg:
            continue
        with open(svg) as f_svg:
            thumbnail = shared.generate_thumbnail(f_svg.read(), svg)
            thumbnail_filename = svg.replace(".colored.", ".thumbnail.")
            with open(thumbnail_filename, "w") as f_thumbnail:
                f_thumbnail.write(thumbnail)
    results_path = Path(results_folder)

    # Move .thumbnail.svg files
    for file in results_path.glob("*.thumbnail.svg"):
        shutil.copy(str(file), folders["thumbnail"])
        file.unlink()

    # Move .colored.svg files
    for file in results_path.glob("*.colored.svg"):
        shutil.copy(str(file), folders["svg"])
        file.unlink()

    # Move .enriched.svg files
    for file in results_path.glob("*.enriched.svg"):
        shutil.copy(str(file), folders["svg"])
        file.unlink()

    # Move .fasta files
    for file in results_path.glob("*.fasta"):
        shutil.copy(str(file), folders["fasta"])
        file.unlink()

    # Move .json files
    for file in results_path.glob("*.json"):
        shutil.copy(str(file), folders["json"])
        file.unlink()


@cli.group("gtrnadb")
def gtrnadb_group():
    """
    Use tRNA templates for structure visualisation.
    """


@gtrnadb_group.command("setup")
def gtrnadb_setup():
    """
    This will copy all the CM files into place so that drawing will not modify
    the data directory.
    """
    rprint(shared.get_r2dt_version_header())
    gtrnadb.setup()


@gtrnadb_group.command("draw")
@click.option(
    "--domain",
    default=False,
    type=click.STRING,
    help="Domain (A for Archaea, B for Bacteria, or E for Eukaryotes)",
)
@click.option(
    "--isotype", default=False, type=click.STRING, help="tRNA isotype, for example Thr"
)
@click.option(
    "--constraint", default=False, is_flag=True, help="Fold insertions using RNAfold"
)
@click.option("--exclusion", default=None)
@click.option("--fold_type", default=None)
@click.option("--quiet", "-q", is_flag=True, default=False)
@click.argument("fasta-input", type=click.Path())
@click.argument("output-folder", type=click.Path())
def gtrnadb_draw(
    fasta_input,
    output_folder,
    domain="",
    isotype="",
    constraint=None,
    exclusion=None,
    fold_type=None,
    quiet=False,
):
    """
    Visualise sequences using GtRNAdb templates.
    """
    # pylint: disable=too-many-arguments,too-many-positional-arguments
    if not quiet:
        rprint(shared.get_r2dt_version_header())
    os.makedirs(output_folder, exist_ok=True)

    fasta_input = shared.sanitise_fasta(fasta_input)

    if domain and isotype:
        core.visualise_trna(
            domain.upper(),
            isotype.capitalize(),
            fasta_input,
            output_folder,
            constraint,
            exclusion,
            fold_type,
            quiet,
        )
    else:
        for trna in gtrnadb.classify_trna_sequences(fasta_input, output_folder):
            core.visualise(
                "gtrnadb",
                fasta_input,
                output_folder,
                trna["id"],
                None,
                constraint,
                exclusion,
                fold_type,
                trna["domain"],
                trna["isotype"],
                trna["start"],
                trna["end"],
                quiet,
            )


@cli.group("rnasep")
def rnasep_group():
    """
    Use RNAse P templates for structure visualisation.
    """


@rnasep_group.command("draw")
@click.option(
    "--constraint", default=False, is_flag=True, help="Fold insertions using RNAfold"
)
@click.option("--exclusion", default=None)
@click.option("--fold_type", default=None)
@click.option(
    "--skip_ribovore_filters",
    default=False,
    is_flag=True,
    help="Ignore ribovore QC checks",
)
@click.option("--quiet", "-q", is_flag=True, default=False)
@click.argument("fasta-input", type=click.Path())
@click.argument("output-folder", type=click.Path())
def rnasep_draw(
    fasta_input,
    output_folder,
    constraint,
    exclusion,
    fold_type,
    quiet,
    skip_ribovore_filters,
):
    """Draw 2D diagrams using RNAse P templates."""
    # pylint: disable=too-many-arguments,too-many-positional-arguments
    if not quiet:
        rprint(shared.get_r2dt_version_header())
    os.makedirs(output_folder, exist_ok=True)

    fasta_input = shared.sanitise_fasta(fasta_input)

    with open(
        shared.get_ribotyper_output(
            fasta_input, output_folder, config.RNASEP_CM_LIBRARY, skip_ribovore_filters
        ),
    ) as f_ribotyper:
        for line in f_ribotyper.readlines():
            rnacentral_id, model_id, _ = line.split("\t")
            core.visualise(
                "rnasep",
                fasta_input,
                output_folder,
                rnacentral_id,
                model_id,
                constraint,
                exclusion,
                fold_type,
                domain=None,
                isotype=None,
                start=None,
                end=None,
                quiet=quiet,
            )


@cli.group("tmrna")
def tmrna_group():
    """
    Use tmRNA templates for structure visualisation.
    """


@tmrna_group.command("draw")
@click.option(
    "--constraint", default=False, is_flag=True, help="Fold insertions using RNAfold"
)
@click.option("--exclusion", default=None)
@click.option("--fold_type", default=None)
@click.option(
    "--skip_ribovore_filters",
    default=False,
    is_flag=True,
    help="Ignore ribovore QC checks",
)
@click.option("--quiet", "-q", is_flag=True, default=False)
@click.argument("fasta-input", type=click.Path())
@click.argument("output-folder", type=click.Path())
def tmrna_draw(
    fasta_input,
    output_folder,
    constraint,
    exclusion,
    fold_type,
    quiet,
    skip_ribovore_filters,
):
    """Draw 2D diagrams using tmRNA templates."""
    # pylint: disable=too-many-arguments,too-many-positional-arguments
    if not quiet:
        rprint(shared.get_r2dt_version_header())
    os.makedirs(output_folder, exist_ok=True)
    with open(
        shared.get_ribotyper_output(
            fasta_input, output_folder, config.TMRNA_CM_LIBRARY, skip_ribovore_filters
        ),
    ) as f_ribotyper:
        for line in f_ribotyper.readlines():
            rnacentral_id, model_id, _ = line.split("\t")
            core.visualise(
                "tmrna",
                fasta_input,
                output_folder,
                rnacentral_id,
                model_id,
                constraint,
                exclusion,
                fold_type,
                domain=None,
                isotype=None,
                start=None,
                end=None,
                quiet=quiet,
            )


@cli.group("crw")
def crw_group():
    """
    Use CRW templates for structure visualisation.
    """


@crw_group.command("draw")
@click.option(
    "--constraint", default=False, is_flag=True, help="Fold insertions using RNAfold"
)
@click.option("--exclusion", default=None)
@click.option("--fold_type", default=None)
@click.option("--quiet", "-q", is_flag=True, default=False)
@click.option(
    "--skip_ribovore_filters",
    default=False,
    is_flag=True,
    help="Ignore ribovore QC checks",
)
@click.argument("fasta-input", type=click.Path())
@click.argument("output-folder", type=click.Path())
def rrna_draw(
    fasta_input,
    output_folder,
    constraint,
    exclusion,
    fold_type,
    quiet,
    skip_ribovore_filters,
):
    """Draw 2D diagrams using CRW templates."""
    # pylint: disable=too-many-arguments,too-many-positional-arguments
    if not quiet:
        rprint(shared.get_r2dt_version_header())
    os.makedirs(output_folder, exist_ok=True)

    fasta_input = shared.sanitise_fasta(fasta_input)

    with open(
        shared.get_ribotyper_output(
            fasta_input, output_folder, config.CRW_CM_LIBRARY, skip_ribovore_filters
        ),
    ) as f_ribotyper:
        for line in f_ribotyper.readlines():
            rnacentral_id, model_id, _ = line.split("\t")
            core.visualise(
                "crw",
                fasta_input,
                output_folder,
                rnacentral_id,
                model_id,
                constraint,
                exclusion,
                fold_type,
                domain=None,
                isotype=None,
                start=None,
                end=None,
                quiet=quiet,
            )


@cli.group("ribovision")
def ribovision_group():
    """
    Use RiboVision templates for structure visualisation.
    """


@ribovision_group.command("draw_lsu")
@click.option(
    "--constraint", default=False, is_flag=True, help="Fold insertions using RNAfold"
)
@click.option("--exclusion", default=None)
@click.option("--fold_type", default=None)
@click.option("--quiet", "-q", is_flag=True, default=False)
@click.option(
    "--skip_ribovore_filters",
    default=False,
    is_flag=True,
    help="Ignore ribovore QC checks",
)
@click.argument("fasta-input", type=click.Path())
@click.argument("output-folder", type=click.Path())
def ribovision_draw_lsu(
    fasta_input,
    output_folder,
    constraint,
    exclusion,
    fold_type,
    quiet,
    skip_ribovore_filters,
):
    """Draw 2D diagrams using LSU templates from RiboVision."""
    # pylint: disable=too-many-arguments,too-many-positional-arguments
    if not quiet:
        rprint(shared.get_r2dt_version_header())
    os.makedirs(output_folder, exist_ok=True)

    fasta_input = shared.sanitise_fasta(fasta_input)

    with open(
        shared.get_ribotyper_output(
            fasta_input,
            output_folder,
            config.RIBOVISION_LSU_CM_LIBRARY,
            skip_ribovore_filters,
        ),
    ) as f_ribotyper:
        for line in f_ribotyper.readlines():
            rnacentral_id, model_id, _ = line.split("\t")
            core.visualise(
                "lsu",
                fasta_input,
                output_folder,
                rnacentral_id,
                model_id,
                constraint,
                exclusion,
                fold_type,
                domain=None,
                isotype=None,
                start=None,
                end=None,
                quiet=quiet,
            )


@ribovision_group.command("draw_ssu")
@click.option(
    "--constraint", default=False, is_flag=True, help="Fold insertions using RNAfold"
)
@click.option("--exclusion", default=None)
@click.option("--fold_type", default=None)
@click.option("--quiet", "-q", is_flag=True, default=False)
@click.option(
    "--skip_ribovore_filters",
    default=False,
    is_flag=True,
    help="Ignore ribovore QC checks",
)
@click.argument("fasta-input", type=click.Path())
@click.argument("output-folder", type=click.Path())
def ribovision_draw_ssu(
    fasta_input,
    output_folder,
    constraint,
    exclusion,
    fold_type,
    quiet,
    skip_ribovore_filters,
):
    """Draw 2D diagrams using SSU templates from RiboVision."""
    # pylint: disable=too-many-arguments,too-many-positional-arguments
    if not quiet:
        rprint(shared.get_r2dt_version_header())
    os.makedirs(output_folder, exist_ok=True)

    fasta_input = shared.sanitise_fasta(fasta_input)

    with open(
        shared.get_ribotyper_output(
            fasta_input,
            output_folder,
            config.RIBOVISION_SSU_CM_LIBRARY,
            skip_ribovore_filters,
        ),
    ) as f_ribotyper:
        for line in f_ribotyper.readlines():
            rnacentral_id, model_id, _ = line.split("\t")
            core.visualise(
                "ssu",
                fasta_input,
                output_folder,
                rnacentral_id,
                model_id,
                constraint,
                exclusion,
                fold_type,
                domain=None,
                isotype=None,
                start=None,
                end=None,
                quiet=quiet,
            )


@cli.group("rfam")
def rfam_group():
    """
    Use Rfam templates for structure visualisation.
    """


@rfam_group.command("blacklisted")
def rfam_blacklist():
    """
    Show all blacklisted families. These include rRNA families as well as
    families that do not have any secondary structure.
    """
    for model in sorted(rfam.blacklisted()):
        rprint(model)


@rfam_group.command("draw")
@click.option(
    "--constraint", default=False, is_flag=True, help="Fold insertions using RNAfold"
)
@click.option("--exclusion", default=None)
@click.option("--fold_type", default=None)
@click.option("--quiet", "-q", is_flag=True, default=False)
@click.option("--rnartist", default=False, is_flag=True)
@click.option("--rscape", default=False, is_flag=True)
@click.argument("rfam_acc", type=click.STRING)
@click.argument("fasta-input", type=click.Path())
@click.argument("output-folder", type=click.Path())
def rfam_draw(
    rfam_acc,
    fasta_input,
    output_folder,
    constraint=None,
    exclusion=None,
    fold_type=None,
    quiet=False,
    rnartist=False,
    rscape=False,
):
    """
    Visualise sequences using the Rfam/R-scape consensus structure as template.

    RFAM_ACCESSION - Rfam family to process (RF00001, RF00002 etc)
    """
    # pylint: disable=too-many-arguments,too-many-positional-arguments
    if not quiet:
        rprint(shared.get_r2dt_version_header())
        rprint(rfam_acc)
    if rnartist and rscape:
        rprint("Please specify only one template type")
        return
    if rnartist:
        template_type = "rnartist"
    elif rscape:
        template_type = "rscape"
    else:
        template_type = "auto"

    fasta_input = shared.sanitise_fasta(fasta_input)

    if rfam.has_structure(rfam_acc):
        rfam.generate_2d(
            rfam_acc,
            output_folder,
            fasta_input,
            constraint,
            exclusion,
            fold_type,
            quiet,
            rfam_template_type=template_type,
        )
    else:
        rprint(f"{rfam_acc} does not have a conserved secondary structure")


@rfam_group.command("validate")
@click.argument("rfam_accession", type=click.STRING)
@click.argument("output", type=click.File("w"))
def rfam_validate(rfam_accession, output):
    """
    Check if the given Rfam accession is one that should be drawn. If so it will
    be output to the given file, otherwise it will not.
    """
    if rfam_accession not in rfam.blacklisted():
        output.write(f"{rfam_accession}\n")


def organise_metadata(output_folder, result_folders):
    """
    Aggregate hits.txt files from all subfolders.
    """
    tsv_folder = os.path.join(output_folder, "results", "tsv")
    os.makedirs(tsv_folder, exist_ok=True)
    with open(os.path.join(tsv_folder, "metadata.tsv"), "w") as f_out:
        for folder in result_folders:
            hits = os.path.join(folder, "hits.txt")
            if not os.path.exists(hits):
                continue
            with open(hits) as f_hits:
                for line in f_hits.readlines():
                    if "gtrnadb" in folder:
                        line = line.replace("PASS", "GtRNAdb")
                    elif "crw" in folder:
                        line = line.replace("PASS", "CRW").replace("FAIL", "CRW")
                    elif "rfam" in folder or "RF00005" in folder:
                        line = line.replace("PASS", "Rfam").replace("FAIL", "Rfam")
                    elif "ribovision-lsu" in folder or "ribovision-ssu" in folder:
                        line = line.replace("PASS", "RiboVision").replace(
                            "FAIL", "RiboVision"
                        )
                    elif "rnasep" in folder:
                        line = line.replace("PASS", "RNAse P Database").replace(
                            "FAIL", "RNAse P Database"
                        )
                    elif "tmrna" in folder:
                        line = line.replace("PASS", "tmRNA Database").replace(
                            "FAIL", "tmRNA Database"
                        )
                    f_out.write(line)


@cli.command()
@click.argument("cm_library", type=click.Path())
def generatemodelinfo(cm_library):
    """
    Helper for generating modelinfo.txt files.
    """
    rprint(shared.get_r2dt_version_header())
    gmi.generate_model_info(cm_library)


def force_draw(
    model_id,
    fasta_input,
    output_folder,
    seq_id,
    constraint=None,
    exclusion=None,
    fold_type=None,
    quiet=False,
):
    """Draw 2D diagrams using a specified template."""
    # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
    if not quiet:
        rprint(shared.get_r2dt_version_header())

    fasta_input = shared.sanitise_fasta(fasta_input)

    model_type = lm.get_model_type(model_id)
    if not model_type:
        rprint("Error: Model not found. Please check model_id")
        return
    if not quiet:
        rprint(
            f"Visualising sequence {seq_id} using the {model_id} model from {model_type}"
        )
    if not os.path.exists(f"{fasta_input}.ssi"):
        runner.run(f"esl-sfetch --index {fasta_input}")

    output = os.path.join(output_folder, model_type.replace("_", "-"))

    if model_type in ["crw", "rfam", "rnasep", "tmrna", "local_data"]:
        core.visualise(
            model_type,
            fasta_input,
            output,
            seq_id,
            model_id,
            constraint,
            exclusion,
            fold_type,
            domain=None,
            isotype=None,
            start=None,
            end=None,
            quiet=quiet,
        )
    elif model_type in ["ribovision_ssu", "ribovision_lsu"]:
        core.visualise(
            model_type.split("_")[1],  # ssu or lsu
            fasta_input,
            output,
            seq_id,
            model_id,
            constraint,
            exclusion,
            fold_type,
            domain=None,
            isotype=None,
            start=None,
            end=None,
            quiet=quiet,
        )
    elif model_type == "gtrnadb":
        domain, isotype = model_id.split("_")
        core.visualise_trna(
            domain,
            isotype,
            fasta_input,
            output,
            constraint,
            exclusion,
            fold_type,
            quiet,
        )
    # organise results into folders
    organise_results(output, output_folder)
    metadata_folder = os.path.join(output_folder, "results", "tsv")
    if not os.path.exists(metadata_folder):
        os.makedirs(metadata_folder)
    label_mapping = {
        "crw": "CRW",
        "gtrnadb": "GtRNAdb",
        "rfam": "Rfam",
        "ribovision_ssu": "RiboVision",
        "ribovision_lsu": "RiboVision",
        "rnasep": "RNAse P database",
        "tmrna": "tmRNA database",
        "local_data": "Local data",
    }
    with open(
        os.path.join(metadata_folder, "metadata.tsv"), "a", encoding="utf-8"
    ) as f_out:
        line = f"{seq_id}\t{model_id}\t{label_mapping[model_type]}\n"
        f_out.write(line)


def _count_overlaps(json_path, threshold=10.0):
    """
    Count nucleotide overlaps in a Traveler JSON file.
    Two nucleotides overlap if their Euclidean distance is less than threshold.
    """
    # pylint: disable=import-outside-toplevel
    import math

    # Check if file exists and is not empty
    if not json_path.exists() or json_path.stat().st_size == 0:
        return float("inf")

    try:
        with open(json_path, "r") as f:
            data = json.load(f)
    except json.JSONDecodeError:
        return float("inf")

    sequence = data["rnaComplexes"][0]["rnaMolecules"][0]["sequence"]
    coords = [(nuc["x"], nuc["y"]) for nuc in sequence]

    overlaps = 0
    n = len(coords)
    for i in range(n):
        for j in range(i + 2, n):  # Skip adjacent nucleotides (i+1)
            dx = coords[i][0] - coords[j][0]
            dy = coords[i][1] - coords[j][1]
            dist = math.sqrt(dx * dx + dy * dy)
            if dist < threshold:
                overlaps += 1
    return overlaps


def _fix_svg_font_size(svg_path, min_font_size=8.0):
    """
    Fix tiny font sizes in SVGs generated from RNArtist templates.
    RNArtist can place nucleotides very close together, causing Traveler
    to use a tiny font. This post-processes the SVG to ensure readable fonts.
    """
    # pylint: disable=import-outside-toplevel,redefined-outer-name,reimported
    import re

    with open(svg_path, "r") as f:
        content = f.read()

    # Find font-size in style attributes (e.g., font-size: 0.828000px)
    def fix_font(match):
        size = float(match.group(1))
        if size < min_font_size:
            return f"font-size: {min_font_size}px"
        return match.group(0)

    content = re.sub(r"font-size:\s*([0-9.]+)px", fix_font, content)

    with open(svg_path, "w") as f:
        f.write(content)


def _templatefree_auto(fasta_input, output_folder, quiet):
    """
    Run R2R, RNApuzzler, and RNArtist, return the layout with fewest overlaps.
    Prefers R2R if overlap counts are equal or very similar.
    """
    # pylint: disable=import-outside-toplevel,too-many-branches,too-many-statements,too-many-locals
    import tempfile

    fasta_input = shared.sanitise_fasta(fasta_input)
    seq_id, sequence, structure = r2r.parse_fasta(fasta_input)

    output_folder = Path(output_folder)
    results_folder = output_folder / "results"

    # Create temp directories for all layouts
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        r2r_output = tmpdir / "r2r"
        rnartist_output = tmpdir / "rnartist"
        puzzler_output = tmpdir / "rnapuzzler"

        # ------------------------------------------------------------------
        # Run R2R
        # ------------------------------------------------------------------
        r2r_folder = r2r_output / "r2r"
        for folder in [r2r_output, r2r_folder]:
            folder.mkdir(exist_ok=True, parents=True)

        r2r.generate_r2r_input_file(sequence, structure, r2r_folder)
        r2r_svg = r2r.run_r2r(r2r_folder)
        rscape_one_line_svg = rfam.convert_rscape_svg_to_one_line(r2r_svg, r2r_folder)
        rfam.convert_rscape_svg_to_traveler(rscape_one_line_svg, r2r_folder)
        scale_coordinates(r2r_folder / "traveler-template.xml", scaling_factor=3)
        r2r.run_traveler(fasta_input, r2r_folder, seq_id)
        organise_results(r2r_folder, r2r_output)
        r2r_json = r2r_output / "results" / "json" / f"{seq_id}.colored.json"

        r2r_overlaps = _count_overlaps(r2r_json)

        if r2r_overlaps == 0:
            # R2R has no overlaps — no need to run the others
            if not quiet:
                rprint("[green]R2R overlaps: 0 -> Using R2R[/green]")
            winner = "R2R"
            winner_output = r2r_output
        else:
            # ----------------------------------------------------------
            # Run RNApuzzler
            # ----------------------------------------------------------
            puzzler_folder_work = puzzler_output / "rnapuzzler"
            for folder in [puzzler_output, puzzler_folder_work]:
                folder.mkdir(exist_ok=True, parents=True)

            try:
                rnapuzzler.run_puzzler_pipeline(
                    fasta_input,
                    puzzler_folder_work,
                    seq_id,
                    sequence,
                    structure,
                )
                organise_results(puzzler_folder_work, puzzler_output)
                puzzler_json = (
                    puzzler_output / "results" / "json" / f"{seq_id}.colored.json"
                )
                puzzler_overlaps = _count_overlaps(puzzler_json)
            except Exception:  # pylint: disable=broad-except
                puzzler_overlaps = float("inf")

            # ----------------------------------------------------------
            # Run RNArtist
            # ----------------------------------------------------------
            rnartist_folder_work = rnartist_output / "rnartist"
            for folder in [rnartist_output, rnartist_folder_work]:
                folder.mkdir(exist_ok=True, parents=True)

            rnartist_obj = RnaArtist(destination=rnartist_folder_work)
            rnartist_obj.fasta_file = fasta_input
            rnartist_obj.seq_label = seq_id
            rnartist_obj.run(rerun=True)

            rnartist_template = rnartist_folder_work / "rnartist-template.xml"
            scale_coordinates(rnartist_template, scaling_factor=3)

            cmd = f"""
            traveler \
                --verbose \
                --target-structure {fasta_input} \
                --template-structure --file-format traveler \
                    {rnartist_folder_work}/rnartist-template.xml {fasta_input} \
                --all {rnartist_folder_work}/{seq_id}
            """
            runner.run(cmd)
            organise_results(rnartist_folder_work, rnartist_output)

            rnartist_svg = rnartist_output / "results" / "svg" / f"{seq_id}.colored.svg"
            if rnartist_svg.exists():
                _fix_svg_font_size(rnartist_svg)

            rnartist_json = (
                rnartist_output / "results" / "json" / f"{seq_id}.colored.json"
            )
            rnartist_overlaps = _count_overlaps(rnartist_json)

            # ----------------------------------------------------------
            # Pick the winner (prefer R2R on ties)
            # ----------------------------------------------------------
            candidates = [
                ("R2R", r2r_overlaps, r2r_output),
                ("RNApuzzler", puzzler_overlaps, puzzler_output),
                ("RNArtist", rnartist_overlaps, rnartist_output),
            ]
            # Sort by overlaps; on ties the original order (R2R first) wins
            candidates.sort(key=lambda c: c[1])
            winner, _, winner_output = candidates[0]

            if not quiet:
                rprint(
                    f"[green]R2R overlaps: {r2r_overlaps}, "
                    f"RNApuzzler overlaps: {puzzler_overlaps}, "
                    f"RNArtist overlaps: {rnartist_overlaps} "
                    f"-> Using {winner}[/green]"
                )

        # Copy winner to output
        output_folder.mkdir(exist_ok=True, parents=True)
        for item in winner_output.iterdir():
            dest = output_folder / item.name
            if item.is_dir():
                if dest.exists():
                    shutil.rmtree(dest)
                shutil.copytree(item, dest)
            else:
                shutil.copy2(item, dest)

        # Write metadata
        tsv_folder = results_folder / "tsv"
        tsv_folder.mkdir(exist_ok=True, parents=True)
        with open(tsv_folder / "metadata.tsv", "w") as f_out:
            f_out.write(f"{seq_id}\t{winner}\t{winner}\n")

        shutil.copyfile(
            fasta_input,
            results_folder / "fasta" / f"{seq_id}.fasta",
        )

        # Remove intermediate RNApuzzler FASTA files from results
        fasta_folder = results_folder / "fasta"
        for name in ("rnapuzzler-template.fasta", "rnapuzzler-target.fasta"):
            (fasta_folder / name).unlink(missing_ok=True)


@cli.command()
@click.argument("fasta-input", type=click.Path())
@click.argument("output-folder", type=click.Path())
@click.option("--quiet", "-q", is_flag=True, default=False)
@click.option("--rnartist", default=False, is_flag=True)
@click.option("--rscape", default=False, is_flag=True)
@click.option(
    "--rnapuzzler",
    "rnapuzzler_flag",
    default=False,
    is_flag=True,
    help="Use RNApuzzler for overlap-free layout (ViennaRNA)",
)
@click.option(
    "--auto/--no-auto",
    default=True,
    help="Run R2R, RNApuzzler, and RNArtist, pick best layout (default)",
)
# pylint: disable=too-many-arguments,too-many-positional-arguments
# pylint: disable=too-many-statements,inconsistent-return-statements,too-many-branches
def templatefree(
    fasta_input, output_folder, rnartist, rscape, rnapuzzler_flag, auto, quiet
):
    """
    Run template-free visualisation using R2R, RNArtist, or RNApuzzler.
    """
    if not quiet:
        rprint(shared.get_r2dt_version_header())

    # Handle auto mode (default unless a specific engine is chosen)
    chosen = sum([rnartist, rscape, rnapuzzler_flag])
    if auto and chosen == 0:
        return _templatefree_auto(fasta_input, output_folder, quiet)

    chosen = sum([rnartist, rscape, rnapuzzler_flag])
    if chosen == 0:
        rscape = True
    if chosen > 1:
        raise ValueError("Please specify only one template type")

    fasta_input = shared.sanitise_fasta(fasta_input)

    if rnartist:
        output_folder = Path(output_folder)
        results_folder = output_folder / "results"
        rnartist_folder = output_folder / "rnartist"

        seq_id, _, _ = r2r.parse_fasta(fasta_input)
        rnartist = RnaArtist(destination=rnartist_folder)
        rnartist.fasta_file = fasta_input
        rnartist.seq_label = seq_id
        rnartist.run(rerun=True)

        # Scale RNArtist coordinates to avoid Traveler numerical issues
        # RNArtist can place nucleotides very close together in tight loops
        rnartist_template = rnartist_folder / "rnartist-template.xml"
        scale_coordinates(rnartist_template, scaling_factor=3)

        cmd = f"""
        traveler \
            --verbose \
            --target-structure {fasta_input} \
            --template-structure --file-format traveler \
                {rnartist_folder}/rnartist-template.xml {fasta_input} \
            --all {rnartist_folder}/{seq_id}
        """
        runner.run(cmd)

        organise_results(rnartist_folder, output_folder)

        # Fix font size in RNArtist SVGs (Traveler uses tiny fonts for tight layouts)
        rnartist_svg = results_folder / "svg" / f"{seq_id}.colored.svg"
        if rnartist_svg.exists():
            _fix_svg_font_size(rnartist_svg)

        tsv_folder = results_folder / "tsv"
        tsv_folder.mkdir(exist_ok=True)
        with open(tsv_folder / "metadata.tsv", "w") as f_out:
            f_out.write(f"{seq_id}\tRNArtist\tRNArtist\n")
        try:
            shutil.copyfile(fasta_input, results_folder / "fasta" / f"{seq_id}.fasta")
            shutil.move(
                rnartist_folder / f"rnartist_{seq_id}.svg",
                results_folder / "svg" / f"{seq_id}.rnartist.svg",
            )
        except FileNotFoundError:
            pass

    if rscape:
        output_folder = Path(output_folder)
        results_folder = output_folder / "results"
        r2r_folder = output_folder / "r2r"
        for folder in [output_folder, results_folder, r2r_folder]:
            folder.mkdir(exist_ok=True, parents=True)

        seq_id, sequence, structure = r2r.parse_fasta(fasta_input)
        r2r.generate_r2r_input_file(sequence, structure, r2r_folder)
        r2r_svg = r2r.run_r2r(r2r_folder)
        rscape_one_line_svg = rfam.convert_rscape_svg_to_one_line(r2r_svg, r2r_folder)
        rfam.convert_rscape_svg_to_traveler(rscape_one_line_svg, r2r_folder)
        scale_coordinates(r2r_folder / "traveler-template.xml", scaling_factor=3)

        has_structure = any(c in structure for c in "()[]{}<>")
        if has_structure:
            r2r.run_traveler(fasta_input, r2r_folder, seq_id)
        else:
            # No base pairs — Traveler segfaults on all-dots structures.
            # Convert the R2R SVG directly into Traveler-compatible format.
            r2r.r2r_svg_to_colored(
                str(r2r_folder / "traveler-template.svg"),
                str(r2r_folder / f"{seq_id}.colored.svg"),
            )

        organise_results(r2r_folder, output_folder)
        tsv_folder = results_folder / "tsv"
        tsv_folder.mkdir(exist_ok=True, parents=True)
        with open(tsv_folder / "metadata.tsv", "w") as f_out:
            f_out.write(f"{seq_id}\tR2R\tR2R\n")
        shutil.copyfile(
            fasta_input,
            results_folder / "fasta" / f"{seq_id}.fasta",
        )

    if rnapuzzler_flag:
        output_folder = Path(output_folder)
        results_folder = output_folder / "results"
        puzzler_folder = output_folder / "rnapuzzler"
        for folder in [output_folder, results_folder, puzzler_folder]:
            folder.mkdir(exist_ok=True, parents=True)

        seq_id, sequence, structure = r2r.parse_fasta(fasta_input)
        rnapuzzler.run_puzzler_pipeline(
            fasta_input, puzzler_folder, seq_id, sequence, structure
        )

        organise_results(puzzler_folder, output_folder)

        # Remove intermediate RNApuzzler FASTA files from results
        fasta_folder = results_folder / "fasta"
        for name in ("rnapuzzler-template.fasta", "rnapuzzler-target.fasta"):
            (fasta_folder / name).unlink(missing_ok=True)

        tsv_folder = results_folder / "tsv"
        tsv_folder.mkdir(exist_ok=True, parents=True)
        with open(tsv_folder / "metadata.tsv", "w") as f_out:
            f_out.write(f"{seq_id}\tRNApuzzler\tRNApuzzler\n")
        shutil.copyfile(
            fasta_input,
            results_folder / "fasta" / f"{seq_id}.fasta",
        )


# =============================================================================
# Stockholm Alignment Processing
# =============================================================================


@cli.command()
@click.argument("stockholm-input", type=click.Path(exists=True))
@click.argument("output-folder", type=click.Path())
@click.option("--quiet", "-q", is_flag=True, default=False)
@click.option(
    "--include-novel",
    is_flag=True,
    default=False,
    help="Also process novelSS_names regions",
)
@click.option(
    "--stitch/--no-stitch",
    default=True,
    help="Automatically stitch output SVGs (default: True)",
)
@click.option(
    "--stitch-output",
    type=click.Path(),
    default=None,
    help="Output path for stitched SVG (default: <output-folder>/stitched.svg)",
)
@click.option(
    "--monochrome/--color",
    default=True,
    help="Monochrome stitched output (default) or preserve colors",
)
@click.option(
    "--color-by",
    type=click.Choice(["none", "structure", "region"], case_sensitive=False),
    default="none",
    help=(
        "Auto-colour panels: 'structure' = unique colour per structureID, "
        "'region' = same colour for all structures sharing a regionID "
        "(overrides --monochrome)"
    ),
)
@click.option(
    "--color-config",
    type=click.Path(exists=True),
    default=None,
    help="TSV file mapping structure/region names to SVG colours (overrides --monochrome)",
)
@click.option(
    "--auto-repair/--no-auto-repair",
    default=False,
    help="Attempt to auto-repair unbalanced bracket structures",
)
@click.option(
    "--all-structures/--named-only",
    default=False,
    help="Include unnamed segments that contain secondary structure (default: named only)",
)
@click.option(
    "--max-unpaired",
    type=int,
    default=30,
    help="Split regions at unpaired stretches longer than N nucleotides (0=disabled)",
)
# pylint: disable=too-many-arguments,too-many-branches,too-many-statements,too-many-locals
# pylint: disable=too-many-positional-arguments
def stockholm(
    stockholm_input,
    output_folder,
    quiet,
    include_novel,
    stitch,
    stitch_output,
    monochrome,
    color_by,
    color_config,
    auto_repair,
    all_structures,
    max_unpaired,
):
    """
    Process a Stockholm alignment with named secondary structure regions.

    Parses the Stockholm file, extracts structural elements and generates
    template-free visualizations for each region. Supports three annotation
    formats:

    **New format** (preferred)::

        #=GC structureID  |..SLI..|..SLII..|...
        #=GC regionID     |..5'UTR.........|..core_protein...

    **Legacy format** (fallback)::

        #=GC knownSS_names  |..SLI..|..SLII..|...

    **Simple alignment** (e.g. Rfam seed)::

        No structureID/regionID required.
        Uses #=GC RF for consensus, #=GF ID/AC for naming.
        Entire alignment treated as one region; stitch auto-skipped.

    Optionally stitches all outputs into a single SVG.

    Example::

        r2dt.py stockholm alignment.stk output_folder
    """
    # pylint: disable=import-outside-toplevel,redefined-outer-name
    from utils import coloring
    from utils import stitch as stitch_module

    if not quiet:
        rprint(shared.get_r2dt_version_header())

    output_folder = Path(output_folder)

    # Process the Stockholm alignment
    processed_regions = stockholm_utils.process_stockholm_alignment(
        stockholm_path=Path(stockholm_input),
        output_folder=output_folder,
        include_novel=include_novel,
        quiet=quiet,
        auto_repair=auto_repair,
        include_all=all_structures,
        max_unpaired=max_unpaired,
        fallback_name=Path(stockholm_input).stem,
    )

    if not processed_regions:
        rprint("[yellow]No regions were successfully processed[/yellow]")
        return

    # Stitch the outputs if requested
    if stitch and len(processed_regions) >= 2:
        if not quiet:
            rprint("\n[blue]Stitching SVG outputs...[/blue]")

        # Collect SVG paths sorted by position
        svg_paths = [Path(r["svg_path"]) for r in processed_regions]

        # Parse SVGs
        panels = []
        for svg_path in svg_paths:
            try:
                panel = stitch_module.parse_svg(svg_path)
                panels.append(panel)
            except (FileNotFoundError, ValueError) as e:
                if not quiet:
                    rprint(f"[yellow]Warning: Could not parse {svg_path}: {e}[/yellow]")

        if len(panels) >= 2:
            # Structure names as per-panel captions.
            # Unnamed regions get empty captions.  When a named region
            # has been split into sub-panels, only the first sub-panel
            # carries the parent's name; the rest get empty captions.
            captions = []
            seen_parents: set[str] = set()
            for r in processed_regions:
                if r.get("is_unnamed"):
                    captions.append("")
                elif r.get("parent_name"):
                    parent = r["parent_name"]
                    if parent not in seen_parents:
                        seen_parents.add(parent)
                        captions.append(parent)
                    else:
                        captions.append("")
                else:
                    captions.append(r["name"])

            # Build region label spans: group consecutive panels by region
            region_labels = []
            if any(r.get("region") for r in processed_regions):
                current_region = None
                span_start = 0
                for idx, r in enumerate(processed_regions):
                    rname = r.get("region") or ""
                    if rname != current_region:
                        if current_region:
                            region_labels.append((current_region, span_start, idx - 1))
                        current_region = rname
                        span_start = idx
                if current_region:
                    region_labels.append(
                        (current_region, span_start, len(processed_regions) - 1)
                    )

            # Build per-panel colours (if requested)
            panel_colors = None
            if color_config:
                panel_colors = coloring.build_color_map(
                    processed_regions, "config", Path(color_config)
                )
            elif color_by in ("structure", "region"):
                panel_colors = coloring.build_color_map(processed_regions, color_by)

            # Stitch
            combined = stitch_module.stitch_svgs(
                panels=panels,
                gap=100,
                glyph_type="break",
                captions=captions,
                monochrome=monochrome and not panel_colors,
                show_gap_labels=True,
                outline=True,
                region_labels=region_labels if region_labels else None,
                panel_colors=panel_colors,
            )

            # Write output
            if stitch_output:
                output_path = Path(stitch_output)
            else:
                output_path = output_folder / "stitched.svg"

            stitch_module.write_svg(combined, output_path)

            if not quiet:
                rprint(f"[green]✓[/green] Stitched SVG written to: {output_path}")

            # Also create an outline-only version for high-level overview
            import copy

            outline_root = copy.deepcopy(combined)
            stitch_module.create_outline_svg(outline_root, stroke_width=4.0)

            outline_path = output_path.with_stem(output_path.stem + "-outline")
            stitch_module.write_svg(outline_root, outline_path)

            if not quiet:
                rprint(f"[green]✓[/green] Outline SVG written to: {outline_path}")

            # Create a clean thumbnail: no annotations, no glyphs, uniform colour
            thumb_panels = []
            for svg_path in svg_paths:
                try:
                    thumb_panels.append(stitch_module.parse_svg(svg_path))
                except (FileNotFoundError, ValueError):
                    pass

            if len(thumb_panels) >= 2:
                thumb_combined = stitch_module.stitch_svgs(
                    panels=thumb_panels,
                    gap=20,
                    glyph_type="none",
                    captions=None,
                    monochrome=not panel_colors,
                    show_gap_labels=False,
                    outline=True,
                    outline_width=4.0,
                    keep_intermediate_labels=False,
                    region_labels=None,
                    panel_colors=panel_colors,
                )
                stitch_module.create_thumbnail_svg(thumb_combined, stroke_width=4.0)

                thumb_path = output_path.with_stem(output_path.stem + "-thumbnail")
                stitch_module.write_svg(thumb_combined, thumb_path)

                if not quiet:
                    rprint(f"[green]✓[/green] Thumbnail SVG written to: {thumb_path}")
        else:
            if not quiet:
                rprint(
                    "[yellow]Not enough valid panels for stitching (need at least 2)[/yellow]"
                )
    elif stitch and len(processed_regions) < 2:
        if not quiet:
            rprint("[yellow]Only one region processed, skipping stitch[/yellow]")

    if not quiet:
        rprint("\n[green]Done![/green]")


# =============================================================================
# SVG Stitching Commands
# =============================================================================


@cli.command()
@click.argument("inputs", nargs=-1, type=click.Path(exists=True))
@click.option(
    "-o", "--output", type=click.Path(), required=True, help="Output SVG file path"
)
@click.option(
    "--gap",
    type=float,
    default=100,
    help="Horizontal gap between panels (default: 100)",
)
@click.option(
    "--sort",
    "sort_by_coords",
    is_flag=True,
    default=False,
    help="Sort by genomic coordinates from filenames",
)
@click.option(
    "--glyph",
    type=click.Choice(["none", "bead", "bar", "break"]),
    default="break",
    help="Join glyph type (default: break)",
)
@click.option(
    "--monochrome/--color",
    default=True,
    help="Monochrome output (default) or preserve colors",
)
@click.option(
    "--captions",
    multiple=True,
    help="Caption for each panel (can specify multiple times)",
)
@click.option(
    "--caption-font-size",
    type=float,
    default=None,
    help="Caption font size (auto-detected if not set)",
)
@click.option(
    "--keep-intermediate-labels",
    is_flag=True,
    default=False,
    help="Keep all 5'/3' labels",
)
@click.option(
    "--no-gap-labels",
    is_flag=True,
    default=False,
    help="Hide nucleotide distance labels",
)
@click.option(
    "--gap-label-font-size",
    type=float,
    default=None,
    help="Font size for gap labels (auto-detected if not set)",
)
@click.option(
    "--no-outline",
    is_flag=True,
    default=False,
    help="Disable connecting outline stroke",
)
@click.option(
    "--outline-color", type=str, default="#cccccc", help="Outline stroke color"
)
@click.option("--outline-width", type=float, default=3.0, help="Outline stroke width")
@click.option(
    "--outline-opacity", type=float, default=0.6, help="Outline opacity (0-1)"
)
@click.option(
    "--anchor-label-font-size",
    type=float,
    default=None,
    help="Font size for 5'/3' labels (auto-detected if not set)",
)
@click.option(
    "--normalize-font-size",
    is_flag=True,
    default=False,
    help="Scale panels to match nucleotide font size of first panel",
)
@click.option("--quiet", "-q", is_flag=True, default=False)
# pylint: disable=too-many-arguments,too-many-locals,too-many-nested-blocks,too-many-positional-arguments
def stitch(
    inputs,
    output,
    gap,
    sort_by_coords,
    glyph,
    monochrome,
    captions,
    caption_font_size,
    keep_intermediate_labels,
    no_gap_labels,
    gap_label_font_size,
    no_outline,
    outline_color,
    outline_width,
    outline_opacity,
    anchor_label_font_size,
    normalize_font_size,
    quiet,
):
    """
    Stitch multiple R2DT SVG diagrams into one combined SVG.

    Arranges panels left-to-right with panel i's 3′ joining panel i+1's 5′.
    """
    # pylint: disable=import-outside-toplevel
    from utils import stitch as stitch_module

    if not quiet:
        rprint(shared.get_r2dt_version_header())

    if len(inputs) < 2:
        rprint("[red]Error: At least 2 input SVG files required[/red]")
        return

    # Parse all input SVGs
    panels = []
    for filepath in inputs:
        try:
            panel = stitch_module.parse_svg(Path(filepath))
            panels.append(panel)
            if not quiet:
                rprint(f"Loaded: {filepath}")
                rprint(f"  viewBox: {panel.viewbox}")
                rprint(f"  5′ anchor: ({panel.anchor_5.x:.2f}, {panel.anchor_5.y:.2f})")
                rprint(f"  3′ anchor: ({panel.anchor_3.x:.2f}, {panel.anchor_3.y:.2f})")
                rprint(
                    f"  genomic coords: {panel.genomic_start:,}-{panel.genomic_end:,}"
                )
        except (FileNotFoundError, ValueError) as e:
            rprint(f"[red]Error: {e}[/red]")
            return

    # Sort by genomic position if requested
    if sort_by_coords:
        panels.sort(key=lambda p: p.sort_key)
        if not quiet:
            rprint("\nSorted order:")
            for i, p in enumerate(panels):
                rprint(
                    f"  {i+1}. {p.filepath.name} ({p.genomic_start:,}-{p.genomic_end:,})"
                )

    # Validate captions
    caption_list = list(captions) if captions else None
    if caption_list and len(caption_list) != len(panels):
        rprint(
            f"[red]Error: --captions requires exactly {len(panels)} values "
            f"(got {len(caption_list)})[/red]"
        )
        return

    # Stitch panels
    result = stitch_module.stitch_svgs(
        panels=panels,
        gap=gap,
        glyph_type=glyph,
        captions=caption_list,
        caption_font_size=caption_font_size,
        caption_pad=8,
        keep_intermediate_labels=keep_intermediate_labels,
        show_gap_labels=not no_gap_labels,
        gap_label_font_size=gap_label_font_size,
        monochrome=monochrome,
        outline=not no_outline,
        outline_color=outline_color,
        outline_width=outline_width,
        outline_opacity=outline_opacity,
        anchor_label_font_size=anchor_label_font_size,
        normalize_font_size=normalize_font_size,
    )

    # Write output
    stitch_module.write_svg(result, Path(output))

    if not quiet:
        rprint(f"\n[green]Created: {output}[/green]")
        rprint(f"  Panels: {len(panels)}")
        rprint(f"  Joins: {len(panels) - 1}")
        if caption_list:
            rprint(f"  Captions: {caption_list}")


# =============================================================================
# Viral Genome Annotation
# =============================================================================


@cli.command("viral-annotate")
@click.argument("genome-fasta", type=click.Path(exists=True))
@click.argument("output-folder", type=click.Path())
@click.option(
    "--stitch-output",
    type=click.Path(),
    default=None,
    help="Path for stitched SVG output",
)
@click.option(
    "--cm-library",
    type=click.Path(),
    default=None,
    help="Path to CM library (default: data/rfam/cms/all.cm)",
)
@click.option(
    "--clanin", type=click.Path(), default=None, help="Path to Rfam.clanin file"
)
@click.option("--cpu", type=int, default=4, help="Number of CPUs for cmscan")
@click.option(
    "--evalue",
    "-E",
    type=float,
    default=None,
    help="E-value threshold (default: use Rfam GA thresholds)",
)
@click.option(
    "--monochrome/--color",
    default=True,
    help="Monochrome output (default) or preserve colors",
)
@click.option("--quiet", "-q", is_flag=True, default=False)
# pylint: disable=too-many-arguments,too-many-branches,too-many-statements,too-many-locals
# pylint: disable=too-many-positional-arguments,too-many-nested-blocks
def viral_annotate(
    genome_fasta,
    output_folder,
    stitch_output,
    cm_library,
    clanin,
    cpu,
    evalue,
    monochrome,
    quiet,
):
    """
    Annotate a viral genome with Rfam ncRNA families and generate diagrams.

    Runs cmscan with Rfam GA thresholds, draws each hit, and stitches results.

    Example:
        r2dt.py viral-annotate genome.fasta output/ --stitch-output stitched.svg
    """
    # pylint: disable=import-outside-toplevel,redefined-outer-name
    from utils import stitch as stitch_module

    if not quiet:
        rprint(shared.get_r2dt_version_header())

    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    # Set default CM library path
    if cm_library is None:
        cm_library = os.path.join(config.RFAM_CM_LIBRARY, "all.cm")

    if not os.path.exists(cm_library):
        rprint(f"[red]Error: CM library not found: {cm_library}[/red]")
        return

    # Check if CM library is indexed (cmpress creates .i1f, .i1i, .i1m, .i1p files)
    cm_i1f = cm_library + ".i1f"
    if not os.path.exists(cm_i1f):
        if not quiet:
            rprint("Indexing CM library with cmpress...")
        runner.run(f"cmpress -F {cm_library} 2>/dev/null")

    # Step 1: Calculate database size
    if not quiet:
        rprint("\n[bold]Step 1: Calculating genome size[/bold]")

    # Try esl-seqstat first, fall back to Python parsing
    total_residues = 0
    seqstat_result = subprocess.run(
        f"esl-seqstat {genome_fasta}",
        shell=True,
        capture_output=True,
        text=True,
        check=False,  # Don't raise on non-zero exit
    )
    if seqstat_result.returncode == 0:
        for line in seqstat_result.stdout.split("\n"):
            # Match both "Total # residues:" and "Total # of residues:"
            if "Total #" in line and "residues:" in line:
                total_residues = int(line.split(":")[1].strip())
                break

    # Fallback: parse FASTA directly with Python
    if total_residues == 0:
        if not quiet:
            rprint("  Using Python FASTA parsing (esl-seqstat not available)")
        with open(genome_fasta) as f:
            for line in f:
                if not line.startswith(">"):
                    total_residues += len(line.strip())

    if total_residues == 0:
        rprint("[red]Error: Could not determine genome size[/red]")
        return

    # Database size in Mb (both strands)
    db_size = (total_residues * 2) / 1_000_000
    if not quiet:
        rprint(f"  Genome size: {total_residues:,} nt")
        rprint(f"  Database size (Z): {db_size:.6f} Mb")

    # Step 2: Run cmscan with Rfam recommended settings
    if not quiet:
        rprint("\n[bold]Step 2: Running cmscan with Rfam GA thresholds[/bold]")

    # Find cmscan executable (may need full path in some environments)
    cmscan_bin = shutil.which("cmscan")
    if cmscan_bin is None:
        # Try common locations
        for path in [
            "/usr/bin/cmscan",
            "/usr/local/bin/cmscan",
            "/opt/infernal/bin/cmscan",
        ]:
            if os.path.exists(path):
                cmscan_bin = path
                break
    if cmscan_bin is None:
        rprint("[red]Error: cmscan not found. Please install Infernal.[/red]")
        return

    tblout_file = output_folder / "cmscan.tblout"
    cmscan_output = output_folder / "cmscan.out"

    # Build cmscan command following Rfam recommendations
    if evalue is not None:
        # Use E-value threshold instead of GA
        threshold_opts = f"-E {evalue}"
        if not quiet:
            rprint(f"  Using E-value threshold: {evalue}")
    else:
        threshold_opts = "--cut_ga"

    cmscan_cmd = (
        f"{cmscan_bin} -Z {db_size:.6f} {threshold_opts} --rfam --nohmmonly "
        f"--tblout {tblout_file} --fmt 2 "
        f"--cpu {cpu} "
    )

    # Add clanin if provided
    if clanin and os.path.exists(clanin):
        cmscan_cmd += f"--clanin {clanin} "

    cmscan_cmd += f"{cm_library} {genome_fasta} > {cmscan_output}"

    if not quiet:
        rprint("  Running: cmscan... (this may take several minutes)")

    exit_code = runner.run(cmscan_cmd)

    if exit_code != 0 or not tblout_file.exists():
        rprint(
            "[red]Error: cmscan failed. Check output/coronavirus/cmscan.out for details.[/red]"
        )
        return

    # Step 3: Parse tblout and filter overlapping hits
    if not quiet:
        rprint("\n[bold]Step 3: Parsing cmscan results[/bold]")

    hits = []
    with open(tblout_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 18:
                continue

            # Format 2 columns: idx, target, accession, query, accession, clan,
            # mdl, mdl_from, mdl_to, seq_from, seq_to, strand, trunc, pass, gc,
            # bias, score, E-value, inc, olp, ...

            try:
                target_name = parts[1]  # Rfam family name
                target_acc = parts[2]  # Rfam accession (e.g., RF00507)
                query_name = parts[3]  # Sequence name
                seq_from = int(parts[9])
                seq_to = int(parts[10])
                strand = parts[11]
                score = float(parts[16])
                evalue = parts[17]
                olp = parts[19] if len(parts) > 19 else "*"

                # Filter: keep only non-overlapped or best-overlapped hits
                # olp = "*" (no overlap), "^" (has overlap but is best), "=" (has better overlap)
                if olp == "=":
                    continue  # Skip hits with better overlapping hits

                # Normalize coordinates (start < end)
                start = min(seq_from, seq_to)
                end = max(seq_from, seq_to)

                hits.append(
                    {
                        "family": target_name,
                        "accession": target_acc,
                        "query": query_name,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "score": score,
                        "evalue": evalue,
                    }
                )

            except (ValueError, IndexError):
                continue

    if not hits:
        rprint("[yellow]No RNA families found in the genome[/yellow]")
        return

    # Sort hits by genomic position
    hits.sort(key=lambda h: h["start"])

    if not quiet:
        rprint(f"  Found {len(hits)} RNA family hits:")
        for h in hits:
            rprint(
                f"    {h['accession']} ({h['family']}): "
                f"{h['start']:,}-{h['end']:,} {h['strand']} score={h['score']:.1f}"
            )

    # Step 4: Extract hit regions and run rfam draw
    if not quiet:
        rprint("\n[bold]Step 4: Generating 2D diagrams for each hit[/bold]")

    # Read genome sequence
    genome_seq = ""
    genome_id = ""
    with open(genome_fasta) as f:
        for line in f:
            if line.startswith(">"):
                genome_id = line[1:].split()[0]
            else:
                genome_seq += line.strip()

    svg_files = []
    rfam_output = output_folder / "rfam"
    rfam_output.mkdir(exist_ok=True)

    for hit in hits:
        # Extract the hit region sequence
        hit_seq = genome_seq[hit["start"] - 1 : hit["end"]]

        if hit["strand"] == "-":
            # Reverse complement
            complement = {
                "A": "U",
                "T": "A",
                "U": "A",
                "G": "C",
                "C": "G",
                "a": "u",
                "t": "a",
                "u": "a",
                "g": "c",
                "c": "g",
            }
            hit_seq = "".join(complement.get(b, b) for b in reversed(hit_seq))

        # Create a FASTA file for this hit
        hit_id = f"{genome_id}_{hit['start']}-{hit['end']}_{hit['strand']}-{hit['accession']}"
        hit_fasta = rfam_output / f"{hit_id}.fasta"

        with open(hit_fasta, "w") as f:
            f.write(f">{hit_id}\n{hit_seq}\n")

        # Run R2DT visualization using the rfam draw command
        try:
            core.visualise(
                "rfam",
                str(hit_fasta),
                str(rfam_output),
                hit_id,
                hit["accession"],
                constraint=False,
                exclusion=None,
                fold_type=None,
                domain=None,
                isotype=None,
                start=None,
                end=None,
                quiet=True,
            )

            # Look for the generated SVG
            svg_pattern = rfam_output / f"{hit_id}.colored.svg"
            if svg_pattern.exists():
                svg_files.append(svg_pattern)
                if not quiet:
                    rprint(f"  [green]✓[/green] {hit['accession']} ({hit['family']})")
            else:
                # Try alternate naming patterns
                for svg in rfam_output.glob(f"*{hit_id}*.svg"):
                    if "colored" in str(svg):
                        svg_files.append(svg)
                        if not quiet:
                            rprint(
                                f"  [green]✓[/green] {hit['accession']} ({hit['family']})"
                            )
                        break
                else:
                    if not quiet:
                        rprint(
                            f"  [yellow]![/yellow] {hit['accession']} - SVG not found"
                        )

        except (FileNotFoundError, ValueError, ET.ParseError) as e:
            if not quiet:
                rprint(f"  [red]✗[/red] {hit['accession']} ({hit['family']}): {e}")

    # Step 5: Stitch SVGs together
    if svg_files and stitch_output:
        if not quiet:
            rprint(f"\n[bold]Step 5: Stitching {len(svg_files)} diagrams[/bold]")

        # Parse and stitch
        panels = []
        for svg_file in svg_files:
            try:
                panel = stitch_module.parse_svg(svg_file)
                panels.append(panel)
            except (FileNotFoundError, ValueError) as e:
                if not quiet:
                    rprint(
                        f"  [yellow]Warning: Could not parse {svg_file}: {e}[/yellow]"
                    )

        if len(panels) >= 2:
            # Sort by genomic position
            panels.sort(key=lambda p: p.sort_key)

            # Generate captions from hit info
            caption_list = [
                f"{h['accession']} ({h['family']})"
                for h in hits
                if any(h["accession"] in str(p.filepath) for p in panels)
            ]

            result = stitch_module.stitch_svgs(
                panels=panels,
                gap=100,
                glyph_type="break",
                captions=caption_list if len(caption_list) == len(panels) else None,
                caption_font_size=36,
                caption_pad=8,
                keep_intermediate_labels=False,
                show_gap_labels=True,
                gap_label_font_size=36,
                monochrome=monochrome,
                outline=True,
                outline_color="#cccccc",
                outline_width=3.0,
                outline_opacity=0.6,
                anchor_label_font_size=48,
                normalize_font_size=True,
            )

            stitch_module.write_svg(result, Path(stitch_output))

            if not quiet:
                rprint(f"\n[green]Created stitched output: {stitch_output}[/green]")
        elif len(panels) == 1:
            # Just copy the single SVG
            shutil.copy(svg_files[0], stitch_output)
            if not quiet:
                rprint(f"\n[green]Copied single diagram to: {stitch_output}[/green]")
        else:
            if not quiet:
                rprint("[yellow]No SVGs available for stitching[/yellow]")

    # Summary
    if not quiet:
        rprint("\n[bold]Summary[/bold]")
        rprint(f"  Genome: {genome_id} ({total_residues:,} nt)")
        rprint(f"  RNA families found: {len(hits)}")
        rprint(f"  Diagrams generated: {len(svg_files)}")
        rprint(f"  Output folder: {output_folder}")


@cli.command()
def list_models():
    """
    List all installed templates.
    """
    rprint(shared.get_r2dt_version_header())
    data = lm.list_models()
    for item in data:
        rprint(item["description"])
    lm.check_unique_descriptions(data)
    with open(os.path.join(config.DATA, "models.json"), "w") as models_file:
        json.dump(data, models_file)


@cli.command()
@click.argument("test_name", required=False, default=None, type=click.STRING)
def test(test_name):
    """
    Run all tests or a special test if provided.
    """
    os.environ["R2DT_KEEP_TEST_RESULTS"] = "1"

    # Discover and run the tests
    loader = unittest.TestLoader()

    if test_name:
        suite = loader.loadTestsFromName(f"tests.tests.{test_name}")
    else:
        suite = loader.discover(".")

    test_runner = unittest.TextTestRunner()
    test_runner.run(suite)


@cli.command()
@click.argument("test_name", required=True, type=click.STRING)
def update_test_examples(test_name):
    """Update test examples for a given test."""
    try:
        class_ = getattr(tests, test_name)
    except AttributeError:
        rprint(f"Error: {test_name} is not found in tests.py")
        return
    test_instance = class_()
    for example_file in test_instance.files:
        rprint(example_file)
        old_filename = os.path.join(
            test_instance.test_results,
            test_instance.test_results_subfolder,
            example_file,
        )
        new_filename = os.path.join(test_instance.precomputed_results, example_file)
        shutil.copyfile(old_filename, new_filename)


@cli.command()
def generatecm():
    """
    Helper for generating covariance models.
    """
    # pylint: disable=redefined-outer-name
    rprint(shared.get_r2dt_version_header())
    for bpseq in glob.glob(f"{config.LOCAL_DATA}/*.bpseq"):
        gcl.convert_bpseq_to_fasta(bpseq)
    for fasta in glob.glob(f"{config.LOCAL_DATA}/*.fasta"):
        rprint(os.path.basename(fasta).replace(".fasta", ""))
        # fasta_no_knots = break_pseudoknots(fasta)
        stockholm = gcl.convert_fasta_to_stockholm(fasta)
        gcl.build_cm(stockholm, config.LOCAL_DATA)
    rprint("Done")


@cli.command()
@click.argument("json_file", type=click.Path())
@click.option("--quiet", "-q", default=False, is_flag=True)
def generate_template(json_file, quiet):
    """Generate an R2DT template from an RNA 2D JSON Schema file."""
    if not quiet:
        rprint(shared.get_r2dt_version_header())
    template = r2djs.SchemaToTemplate(json_file)
    if not quiet:
        rprint(f"Created a new {template}")


# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
# pylint: disable=too-many-branches,too-many-statements
def _run_multichain_pdb(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    ctx,
    file_path,
    actual_format,
    structure_id,
    output_path,
    chains,
    rnapuzzler_flag,
    simulate_model,
    simulate_seed,
    compare,
    quiet,
    model_file=None,
    model_chains=None,
    _skip_partition=False,
):
    """Build a single combined 2D diagram from multiple RNA chains.

    Extracts the selected chains, concatenates them into one label space
    (auto-ordering unless an explicit order is given), lays out the nested pairs
    via the templatefree pipeline, then post-processes the SVG to break the
    phantom backbone bond at each chain junction, label per-chain 5'/3' termini,
    and draw crossing inter-chain pairs as an overlay.  The concatenation
    metadata is also written to ``<id>.multichain.json`` for the
    model-comparison draw.

    In ``--compare`` mode with an explicit ``--chains``, the requested chains
    may not be the whole story: the reference structure could have other RNA
    chains that interact with them (widen the display to the whole group —
    see ``multichain.partition_components``) or other chains entirely
    unrelated to them (offer a chain picker instead of silently dropping
    them). ``_skip_partition`` is for the sibling reference-only pages this
    generates for the picker — they already name an exact chain set and
    shouldn't recurse into building their own siblings.
    """
    # pylint: disable=import-outside-toplevel
    from utils import multichain

    if actual_format != "cif":
        rprint(
            "[red]Multi-chain mode requires an mmCIF structure "
            "(use --format cif or provide a .cif file).[/red]"
        )
        return

    if chains:
        chain_list = [c.strip() for c in chains.split(",") if c.strip()]
        auto_order = False  # explicit order provided
    else:  # --all-chains
        chain_list = None
        auto_order = True

    extraction_dir = output_path / "extraction"

    score_chains = None
    chain_views = None
    if compare and chain_list and not _skip_partition:
        components = multichain.partition_components(
            str(file_path), str(extraction_dir), quiet=True
        )
        if components:
            requested = set(chain_list)
            default_component = next(
                (c for c in components if requested & set(c)), None
            )
            if default_component:
                if set(default_component) != requested:
                    # The requested chains interact with others — widen the
                    # display to the whole group; scoring still only covers
                    # what was actually requested. This can happen even with
                    # only one component total (e.g. a 2-chain dimer where
                    # both chains interact and only one was requested).
                    score_chains = chain_list
                    chain_list = default_component
                    auto_order = True
                # A chain picker is only needed if there's something else in
                # the structure *not* part of the default view.
                other_components = [c for c in components if c != default_component]
                if other_components:
                    chain_views = _build_chain_views(
                        ctx,
                        file_path,
                        actual_format,
                        structure_id,
                        output_path,
                        default_component,
                        other_components,
                        rnapuzzler_flag,
                        quiet,
                    )

    result = multichain.assemble(
        str(file_path),
        str(extraction_dir),
        chains=chain_list,
        auto_order=auto_order,
        quiet=quiet,
    )
    if result is None:
        rprint("[red]Error: multi-chain extraction failed[/red]")
        return

    counts = result.counts()
    if not quiet:
        rprint(
            f"Chains: {'+'.join(result.order)}  "
            f"({counts['length']} nt, {counts['nested']} nested pairs, "
            f"{counts['crossing']} crossing)"
        )

    # Record the concatenation metadata for later phases.
    meta_path = output_path / f"{structure_id}.multichain.json"
    with open(meta_path, "w") as meta_file:
        json.dump(
            {
                "structure_id": structure_id,
                "order": result.order,
                "boundaries": [
                    {"chain": c, "start": s, "end": e} for c, s, e in result.boundaries
                ],
                "sequence": result.sequence,
                "dot_bracket": result.dot_bracket,
                "nested_pairs": result.nested_pairs,
                "crossing_pairs": result.crossing_pairs,
                "counts": counts,
            },
            meta_file,
            indent=2,
        )

    # Write the combined FASTA (sequence + nested dot-bracket) and lay it out.
    fasta_path = output_path / f"{structure_id}.fasta"
    with open(fasta_path, "w") as fasta_file:
        fasta_file.write(f">{structure_id}\n{result.sequence}\n{result.dot_bracket}\n")

    results_folder = output_path / "results"
    if not quiet:
        rprint("Generating combined 2D diagram (templatefree)...")
    # Use auto mode (R2R + RNApuzzler + RNArtist, keep the fewest-overlap layout)
    # unless the caller forced RNApuzzler. Passing rscape=True would set chosen==1
    # and bypass the auto overlap-minimising selection, forcing R2R (which is not
    # overlap-free) -- so leave rscape False here.
    ctx.invoke(
        templatefree,
        fasta_input=str(fasta_path),
        output_folder=str(results_folder),
        rnartist=False,
        rscape=False,
        rnapuzzler_flag=rnapuzzler_flag,
        quiet=quiet,
    )

    svg_dir = results_folder / "results" / "svg"
    svg_path = next(
        (
            p
            for p in (
                svg_dir / f"{structure_id}.colored.svg",
                svg_dir / f"{structure_id}.svg",
            )
            if p.exists()
        ),
        None,
    )
    if svg_path is not None:
        # Break phantom junctions, add per-chain 5'/3' termini, draw the
        # crossing-pair overlay.
        multichain.postprocess_combined_svg(
            str(svg_path),
            result.boundaries,
            result.nested_pairs,
            result.crossing_pairs,
            quiet=quiet,
        )
        if not quiet:
            rprint(f"[green]Success! Combined SVG created: {svg_path}[/green]")

        if simulate_model and not compare:
            _emit_simulated_model_panel(
                svg_path, structure_id, result, simulate_seed, quiet
            )
        if compare:
            _emit_compare_viewer(
                ctx,
                file_path,
                actual_format,
                structure_id,
                output_path,
                result,
                simulate_seed,
                quiet,
                model_file=model_file,
                model_chains=model_chains,
                score_chains=score_chains,
                chain_views=chain_views,
                rnapuzzler_flag=rnapuzzler_flag,
            )
    elif not quiet:
        rprint("[yellow]Diagram generation completed. Check output folder.[/yellow]")


def _chain_group_slug(chain_ids):
    return re.sub(r"[^A-Za-z0-9]+", "-", "-".join(chain_ids)).strip("-") or "chain"


# pylint: disable=too-many-arguments,too-many-positional-arguments
def _build_chain_views(
    ctx,
    file_path,
    actual_format,
    structure_id,
    output_path,
    default_component,
    other_components,
    rnapuzzler_flag,
    quiet,
):
    """Generate a reference-only sibling page for each of ``other_components``
    — RNA chains of the same structure that don't interact with the ones
    shown by default — and return the ``chainViews`` list (the default page
    plus every sibling) for the reference panel's chain picker.

    Each sibling is a real (non-simulated) self-compare: the structure
    against itself, restricted to that component's chains. This reuses the
    existing compare-viewer machinery to get a real interactive 2D+3D page
    with no scoring semantics attached — INF=1.000/RMSD=0.00 there are the
    expected, meaningless artifact of comparing a structure to itself, not a
    real prediction score.
    """
    views = [
        {
            "label": f"Chain {'+'.join(default_component)} (current)",
            "url": "./index.html",
            "current": True,
        }
    ]
    for comp in other_components:
        slug = _chain_group_slug(comp)
        sib_dir = output_path / f"chain-{slug}"
        sib_dir.mkdir(parents=True, exist_ok=True)
        _run_multichain_pdb(
            ctx,
            file_path,
            actual_format,
            structure_id,
            sib_dir,
            ",".join(comp),
            rnapuzzler_flag,
            simulate_model=False,
            simulate_seed=2,
            compare=True,
            quiet=quiet,
            model_file=file_path,
            model_chains=",".join(comp),
            _skip_partition=True,
        )
        views.append(
            {
                # Sibling of the parent output dir, not of viewer/ itself — the
                # current page lives one level down, at <output>/viewer/index.html.
                "label": f"Chain {'+'.join(comp)}",
                "url": f"../chain-{slug}/viewer/index.html",
                "current": False,
            }
        )
    return views


def _emit_simulated_model_panel(svg_path, structure_id, result, seed, quiet):
    """TESTING aid: draw a reference/model base-pair diff on the reference
    layout, using a randomly perturbed copy of the reference pairs as a stand-in
    model (see multichain.simulate_model_pairs)."""
    # pylint: disable=import-outside-toplevel
    from utils import multichain

    ref_pairs = result.nested_pairs + result.crossing_pairs
    model_pairs = multichain.simulate_model_pairs(
        ref_pairs, len(result.sequence), seed=seed
    )
    model_svg = svg_path.with_name(f"{structure_id}.model.svg")
    multichain.render_model_panel(
        str(svg_path), str(model_svg), ref_pairs, model_pairs, quiet=quiet
    )
    if not quiet:
        matched, lost, added = multichain.diff_pairs(ref_pairs, model_pairs)
        rprint(
            f"[cyan]Simulated model panel: {model_svg} "
            f"({len(matched)} matched, {len(lost)} missing, "
            f"{len(added)} model-only)[/cyan]"
        )


def _ensure_cif(path, workdir, quiet):
    """Return an mmCIF path for ``path``, converting a PDB via BioPython.

    ``utils.multichain`` reads structures through FR3D's mmCIF reader, so a
    predicted model supplied as ``.pdb`` is converted first.
    """
    if path.suffix.lower() == ".cif":
        return path
    # pylint: disable=import-outside-toplevel
    import warnings

    from Bio.PDB import MMCIFIO, PDBParser

    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / f"{path.stem}.cif"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = PDBParser(QUIET=True).get_structure(path.stem, str(path))
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(str(out))
    if not quiet:
        rprint(f"[dim]Converted model {path.name} -> {out.name}[/dim]")
    return out


# CASP16 assessment colour scheme (Das et al.): the experimental structure is
# coloured by molecule type — RNA green, DNA orange, ligand red, protein yellow —
# and the best predicted model is overlaid in dark blue.
_CASP_RNA_GREEN = {"r": 0, "g": 166, "b": 81}
_CASP_DNA_ORANGE = {"r": 247, "g": 148, "b": 30}
_CASP_LIGAND_RED = {"r": 237, "g": 28, "b": 36}
_CASP_PROTEIN_YELLOW = {"r": 255, "g": 204, "b": 0}
_CASP_MODEL_BLUE = {"r": 26, "g": 58, "b": 140}


# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
# pylint: disable=too-many-statements,too-many-branches
def _layout_multichain_structure(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    ctx,
    *,
    structure_id: str,
    sequence: str,
    dot_bracket: str,
    boundaries,
    nested_pairs,
    crossing_pairs,
    layout_dir: Path,
    rnapuzzler_flag: bool,
    quiet: bool,
):
    """Run templatefree on a multi-chain FASTA; return (colored_json, colored_svg) or None."""
    # pylint: disable=import-outside-toplevel
    from utils import multichain

    layout_dir = Path(layout_dir)
    layout_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = layout_dir / f"{structure_id}.fasta"
    fasta_path.write_text(
        f">{structure_id}\n{sequence}\n{dot_bracket}\n", encoding="utf-8"
    )
    results_folder = layout_dir / "results"
    ctx.invoke(
        templatefree,
        fasta_input=str(fasta_path),
        output_folder=str(results_folder),
        rnartist=False,
        rscape=False,
        rnapuzzler_flag=rnapuzzler_flag,
        quiet=quiet,
    )
    svg_dir = results_folder / "results" / "svg"
    json_dir = results_folder / "results" / "json"
    svg_path = next(
        (
            p
            for p in (
                svg_dir / f"{structure_id}.colored.svg",
                svg_dir / f"{structure_id}.svg",
            )
            if p.exists()
        ),
        None,
    )
    json_path = json_dir / f"{structure_id}.colored.json"
    if svg_path is None or not json_path.exists():
        return None
    multichain.postprocess_combined_svg(
        str(svg_path),
        boundaries,
        nested_pairs,
        crossing_pairs,
        quiet=quiet,
    )
    return json_path, svg_path


def _emit_compare_viewer(  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-statements,too-many-branches
    ctx,
    source_structure_path,
    actual_format,
    structure_id,
    output_path,
    result,
    simulate_seed,
    quiet,
    model_file=None,
    model_chains=None,
    score_chains=None,
    chain_views=None,
    rnapuzzler_flag=False,
):
    """Assemble the interactive 3-panel compare page.

    Reference 2D and the default model 2D panel share the reference's combined
    layout (approach B). A second model layout (``model-own/``) is also
    generated from the model's own structure so the viewer can switch on
    demand. Each panel's data lives under ``ref/`` / ``model/`` / ``model-own/``
    next to ``index.html``.
    """
    # pylint: disable=import-outside-toplevel
    from utils import multichain

    results_folder = output_path / "results"
    colored_json = results_folder / "results" / "json" / f"{structure_id}.colored.json"
    colored_svg = results_folder / "results" / "svg" / f"{structure_id}.colored.svg"
    if not colored_json.exists():
        if not quiet:
            rprint("[yellow]No colored JSON; cannot build compare viewer.[/yellow]")
        return

    n = len(result.sequence)

    # score_offsets[k] = the absolute position in result's (possibly widened)
    # label space that the model's k-th residue (in its own, narrower, local
    # order) corresponds to. Identity when unwidened, so every computation
    # below that uses it is a no-op in the common case.
    if score_chains:
        by_chain = {cid: (start, end) for cid, start, end in result.boundaries}
        score_offsets = [pos for cid in score_chains for pos in range(*by_chain[cid])]
    else:
        score_offsets = list(range(n))
    score_positions = set(score_offsets)
    score_n = len(score_offsets)

    def _remap_pairs(pairs):
        return [(score_offsets[i], score_offsets[j]) for i, j in pairs]

    def _remap_all_pairs(pairs):
        return [(score_offsets[i], score_offsets[j], fam) for i, j, fam in pairs]

    api_data = viewer_export.build_multichain_api_data(
        colored_json,
        result.chain_of,
        result.auth_of,
        structure_id,
        colored_svg_path=colored_svg if colored_svg.exists() else None,
    )
    # Model panel only: grey out nucleotides outside score_positions (e.g. a
    # widened reference's unmodeled second chain) using the existing
    # unobserved/unresolved-residue dimming path (backbone + text), so it
    # reads as "not part of this model" rather than looking like normal,
    # clickable model content. The reference panel keeps the plain api_data
    # — those residues are real reference data, nothing to greyed there.
    unscored_labels = sorted(p + 1 for p in range(n) if p not in score_positions)
    model_api_data = (
        {**api_data, "unobserved_label_seq_ids": unscored_labels}
        if unscored_labels
        else api_data
    )
    ref_fr3d = viewer_export.build_pairs_fr3d_data(
        result.nested_pairs,
        result.crossing_pairs,
        result.sequence,
        result.auth_of,
        structure_id,
        all_pairs=result.all_pairs,
    )

    # Full pairs (unscoped) for the reference panel's own display/list — shown
    # as-is even when widened, so the diagram honestly reflects the whole
    # structure. Scoped for the *scored* diff below, so a wider display never
    # changes the official matched/lost/added counts.
    ref_pairs = result.nested_pairs + result.crossing_pairs
    scoped_ref_pairs = [
        (i, j) for i, j in ref_pairs if i in score_positions and j in score_positions
    ]

    model_result = None
    model_cif = None
    if model_file:
        # Real predicted model: extract its base pairs and place them on the
        # reference coordinates (approach B).  Requires the same sequence in
        # the same chain order (co-indexed); we validate the length here.
        model_id = Path(model_file).stem
        model_cif = _ensure_cif(Path(model_file), output_path / "downloads", quiet)
        chain_list = (
            [c.strip() for c in model_chains.split(",") if c.strip()]
            if model_chains
            else None
        )
        model_result = multichain.assemble(
            str(model_cif),
            str(output_path / "model_extraction"),
            chains=chain_list,
            auto_order=False,
            quiet=quiet,
        )
        if model_result is None:
            rprint("[red]Error: could not extract base pairs from the model[/red]")
            return
        if len(model_result.sequence) != score_n:
            rprint(
                f"[red]Model/reference length mismatch "
                f"({len(model_result.sequence)} vs {score_n}); the two must share "
                f"the same sequence in the same chain order. Aborting compare.[/red]"
            )
            return
        # Positioned onto result's label space (score_offsets is identity
        # unless the reference display was widened beyond score_chains).
        model_nested = _remap_pairs(model_result.nested_pairs)
        model_crossing = _remap_pairs(model_result.crossing_pairs)
        model_all_pairs = _remap_all_pairs(model_result.all_pairs)
        model_is_simulated = False
    else:
        model_id = f"{structure_id}_model"
        model_pairs = multichain.simulate_model_pairs(ref_pairs, n, seed=simulate_seed)
        model_nested, model_crossing = multichain.max_noncrossing(model_pairs, n)
        # The simulated model only perturbs canonical pairs.
        model_all_pairs = [(i, j, "cWW") for i, j in model_pairs]
        model_is_simulated = True

    # Interaction Network Fidelity of the model against the reference. Scoped
    # to score_chains (a no-op filter when unwidened) so a wider *display*
    # never inflates or deflates the score with pairs the model never
    # attempted (e.g. a second, unpredicted chain of a dimer).
    scoped_ref_all_pairs = [
        (i, j, fam)
        for i, j, fam in result.all_pairs
        if i in score_positions and j in score_positions
    ]
    inf_metrics = multichain.compute_inf(scoped_ref_all_pairs, model_all_pairs)

    model_fr3d = viewer_export.build_pairs_fr3d_data(
        model_nested,
        model_crossing,
        result.sequence,
        result.auth_of,
        model_id,
        all_pairs=model_all_pairs,
    )

    viewer_dir = output_path / "viewer"
    viewer_dir.mkdir(exist_ok=True)
    structure_name = f"{structure_id}.{actual_format}"
    # Keep only the analysed RNA chain(s) in the 3D reference so a deposited
    # entry's extra chains (other RNA copies, ligands, whole proteins) don't
    # render as additional structures beside the one the model overlays.
    if not multichain.write_structure_chains(
        str(source_structure_path),
        list(result.order),
        str(viewer_dir / structure_name),
        quiet=quiet,
    ):
        shutil.copyfile(source_structure_path, viewer_dir / structure_name)
    _copy_viewer_assets(viewer_dir)

    # Superimpose the real predicted model onto the reference so the 3D pane can
    # overlay both structures (approach C: pre-aligned coordinates — the pinned
    # pdbe-molstar bundle exposes no in-browser transform).  Reference and model
    # are co-indexed (approach B), giving an exact per-residue atom correspondence.
    overlays = []
    superpose_rmsd = None
    superpose_n = None
    if model_file and not model_is_simulated:
        # ref_index walks score_offsets (not 0..n-1) so a widened result still
        # only superposes on the chains the model actually corresponds to.
        ref_index = [
            (
                (result.chain_of[p], result.auth_of[p])
                if p < len(result.auth_of) and result.auth_of[p] is not None
                else None
            )
            for p in score_offsets
        ]
        model_index = [
            (
                (model_result.chain_of[i], model_result.auth_of[i])
                if i < len(model_result.auth_of) and model_result.auth_of[i] is not None
                else None
            )
            for i in range(score_n)
        ]
        aligned_name = f"{model_id}.aligned.cif"
        aligned_path = viewer_dir / aligned_name
        sp = multichain.superpose_model_onto_reference(
            str(source_structure_path),
            str(model_cif),
            ref_index,
            model_index,
            str(aligned_path),
            quiet=quiet,
        )
        if sp is None:
            # Too few correspondences to fit — show the model in its own frame.
            shutil.copyfile(str(model_cif), str(aligned_path))
        else:
            superpose_rmsd, superpose_n = sp
        # The model has its own author numbering, so the 3D selection needs a
        # label→(auth, chain) map built from the model (not the reference).
        # Keyed by both shared-layout labels (score_offsets[i]+1) and own-layout
        # labels (i+1) so either model 2D mode can drive the overlay.
        model_label_to_auth = {}
        model_label_to_chain = {}
        for i in range(score_n):
            shared_key = str(score_offsets[i] + 1)
            own_key = str(i + 1)
            if i < len(model_result.auth_of) and model_result.auth_of[i] is not None:
                model_label_to_auth[shared_key] = model_result.auth_of[i]
                model_label_to_auth[own_key] = model_result.auth_of[i]
            if i < len(model_result.chain_of):
                model_label_to_chain[shared_key] = model_result.chain_of[i]
                model_label_to_chain[own_key] = model_result.chain_of[i]
        # label-maps.json lives in model/ alongside that panel's own data (the
        # panels' ref/ and model/ dirs are created below, but this write runs
        # first, so make sure model/ exists here too).
        model_dir = viewer_dir / "model"
        model_dir.mkdir(exist_ok=True)
        (model_dir / "label-maps.json").write_text(
            json.dumps(
                {
                    "labelToAuth": model_label_to_auth,
                    "labelToChain": model_label_to_chain,
                }
            )
        )
        overlays.append(
            {
                "structureId": model_id,
                "structureUrl": aligned_name,
                "structureFormat": "cif",
                # Best predicted model: dark blue (CASP16 scheme).
                "baseColor": _CASP_MODEL_BLUE,
                "baseUrl": "model/",
            }
        )

    chains_label = "+".join(result.order)
    model_tag = "model, simulated" if model_is_simulated else "model"

    # The viewer uses structureId as the 2D plugin's pdbId, which is interpolated
    # into CSS selectors (e.g. querySelector('.rnaTopoSvg_<id>'),
    # '.rnaview_<id>_<n>') by both the plugin and the glue code. A dot in the id
    # (e.g. from a "<model>.trim.pdb" filename → stem "<model>.trim") makes the
    # selector parse the '.' as a *new class* (".rnaTopoSvg_x.trim" = class
    # rnaTopoSvg_x AND class trim), so those lookups match nothing and that panel's
    # selection/click path silently breaks. Sanitise the id to [A-Za-z0-9_].
    def _safe_id(value):
        return re.sub(r"[^A-Za-z0-9_]", "_", str(value))

    # Per-panel base-pair-list annotation: classify each listed pair against the
    # *other* structure's pair set (position keys in the shared label space).
    # Reference list → TP (also in model) / FN (missed); model list → TP / FP.
    def _pair_keys(pairs):
        return sorted({f"{min(i, j) + 1}_{max(i, j) + 1}" for i, j in pairs})

    # Compare over *all* families (cWW + non-canonical), by residue-pair
    # position, so the base-pair list's TP/FP/FN badges are correct for the
    # non-canonical pairs now shown too — not just the cWW backbone.
    ref_pair_keys = _pair_keys((i, j) for i, j, _ in result.all_pairs)
    model_pair_keys = _pair_keys((i, j) for i, j, _ in model_all_pairs)

    # Write each panel's data as sibling JSON files instead of inlining it into
    # index.html: the client already fetches api.json/fr3d.json/lbn.json from a
    # panel's own baseUrl when they aren't passed inline (resolvePanelData() in
    # r2dt-2d-3d-viewer.js), so this is a pure generation-side change. Keeps
    # index.html small and the raw data independently fetchable/cacheable.
    ref_dir = viewer_dir / "ref"
    model_dir = viewer_dir / "model"
    ref_dir.mkdir(exist_ok=True)
    model_dir.mkdir(exist_ok=True)
    ref_lbn = lbn_export.build_lbn_data(api_data, ref_fr3d)
    model_lbn = lbn_export.build_lbn_data(model_api_data, model_fr3d)
    (ref_dir / "api.json").write_text(json.dumps(api_data))
    (ref_dir / "fr3d.json").write_text(json.dumps(ref_fr3d))
    (ref_dir / "lbn.json").write_text(json.dumps(ref_lbn))
    (model_dir / "api.json").write_text(json.dumps(model_api_data))
    (model_dir / "fr3d.json").write_text(json.dumps(model_fr3d))
    (model_dir / "lbn.json").write_text(json.dumps(model_lbn))
    # Also drop a root-level copy of the reference's data (duplicated from
    # ref/): the shared Mol* pane below always links panelIndex 0 (reference)
    # and resolves its own api.json/fr3d.json relative to the viewer root, not
    # the panel's baseUrl, so those two files need to exist there too.
    (viewer_dir / "api.json").write_text(json.dumps(api_data))
    (viewer_dir / "fr3d.json").write_text(json.dumps(ref_fr3d))
    # otherPairKeys is the *other* panel's base-pair position keys, used only
    # to classify this panel's own list as TP/FP/FN -- same size concern as
    # apiData/fr3dData, so it goes in bp-compare.json alongside them rather
    # than inline in the panel dict.
    (ref_dir / "bp-compare.json").write_text(json.dumps(model_pair_keys))
    (model_dir / "bp-compare.json").write_text(json.dumps(ref_pair_keys))

    # Independent model layout (always for real models) so the viewer can
    # switch the right-hand 2D between shared-reference and own coordinates.
    model_own_layout = None
    if model_result is not None and not model_is_simulated:
        if not quiet:
            rprint("Generating model's own 2D layout (templatefree)...")
        own_layout = _layout_multichain_structure(
            ctx,
            structure_id=_safe_id(model_id),
            sequence=model_result.sequence,
            dot_bracket=model_result.dot_bracket,
            boundaries=model_result.boundaries,
            nested_pairs=model_result.nested_pairs,
            crossing_pairs=model_result.crossing_pairs,
            layout_dir=output_path / "model_layout",
            rnapuzzler_flag=rnapuzzler_flag,
            quiet=quiet,
        )
        if own_layout is None:
            if not quiet:
                rprint(
                    "[yellow]Model own-layout generation failed; "
                    "shared layout only.[/yellow]"
                )
        else:
            own_json, own_svg = own_layout
            own_api = viewer_export.build_multichain_api_data(
                own_json,
                model_result.chain_of,
                model_result.auth_of,
                _safe_id(model_id),
                colored_svg_path=own_svg,
            )
            own_fr3d = viewer_export.build_pairs_fr3d_data(
                model_result.nested_pairs,
                model_result.crossing_pairs,
                model_result.sequence,
                model_result.auth_of,
                _safe_id(model_id),
                all_pairs=model_result.all_pairs,
            )
            # Ref pair keys remapped into the model's native 1..score_n space
            # so TP/FP/FN badges stay correct on the own-layout pair list.
            inv = {score_offsets[i]: i for i in range(score_n)}
            ref_keys_own = _pair_keys(
                (inv[i], inv[j])
                for i, j, _ in scoped_ref_all_pairs
                if i in inv and j in inv
            )
            own_dir = viewer_dir / "model-own"
            own_dir.mkdir(exist_ok=True)
            own_lbn = lbn_export.build_lbn_data(own_api, own_fr3d)
            (own_dir / "api.json").write_text(json.dumps(own_api))
            (own_dir / "fr3d.json").write_text(json.dumps(own_fr3d))
            (own_dir / "lbn.json").write_text(json.dumps(own_lbn))
            (own_dir / "bp-compare.json").write_text(json.dumps(ref_keys_own))
            # Dual-key maps (same as model/) so 3D click-through works.
            if (model_dir / "label-maps.json").is_file():
                shutil.copyfile(
                    model_dir / "label-maps.json", own_dir / "label-maps.json"
                )
            to_shared = {str(i + 1): score_offsets[i] + 1 for i in range(score_n)}
            from_shared = {str(score_offsets[i] + 1): i + 1 for i in range(score_n)}
            model_own_layout = {
                "baseUrl": "model-own/",
                "labelBridge": {"toShared": to_shared, "fromShared": from_shared},
            }
            if not quiet:
                rprint(f"[green]Model own layout ready: {own_svg}[/green]")

    panels = [
        {
            "title": f"{structure_id} (reference)",
            "subtitle": f"2D · chains {chains_label}",
            "structureId": _safe_id(structure_id),
            "chainId": "",
            "baseUrl": "ref/",
            "structureUrl": f"../{structure_name}",
            "structureFormat": actual_format,
            "bpCompare": {"role": "reference"},
        },
        {
            "title": f"{model_id} ({model_tag})",
            "subtitle": f"2D · chains {chains_label}",
            "structureId": _safe_id(model_id),
            "chainId": "",
            "baseUrl": "model/",
            "structureUrl": f"../{structure_name}",
            "structureFormat": actual_format,
            "bpCompare": {"role": "model"},
        },
    ]
    if model_own_layout:
        panels[1]["layoutModes"] = {
            "shared": {"baseUrl": "model/", "labelBridge": None},
            "own": model_own_layout,
        }
        panels[1]["defaultLayoutMode"] = "shared"
    if chain_views:
        panels[0]["chainViews"] = chain_views
    molstar = {
        "panelIndex": 0,
        "baseUrl": ".",
        "structureId": structure_id,
        "structureUrl": structure_name,
        "structureFormat": actual_format,
    }
    if overlays:
        # Experimental structure: RNA green (CASP16 scheme). These are RNA targets,
        # so a single RNA colour is used; DNA/ligand/protein colours are defined
        # above for when mixed-molecule targets get per-component colouring.
        molstar["baseColor"] = _CASP_RNA_GREEN
        molstar["overlays"] = overlays
        # Give the 2D panels' click-to-highlight colour the same per-structure
        # colour as the 3D pane (reference green / model blue), so a 2D
        # selection reads as the same colour in both views. Only when a real
        # 3D overlay exists — otherwise the 3D pane has no per-structure
        # colour to match.
        panels[0]["baseColor"] = _CASP_RNA_GREEN
        panels[1]["baseColor"] = _CASP_MODEL_BLUE
    heading = f"{structure_id} — reference vs {model_tag}"
    html_path = viewer_html.render_compare(
        viewer_dir,
        page_title=f"{structure_id} — reference vs model",
        heading=heading,
        subtitle=f"chains {chains_label} · shared 3D",
        panels=panels,
        molstar=molstar,
        metrics=inf_metrics,
    )

    # Structured metrics for the batch dashboard (avoids parsing stdout).
    # Scoped, like INF above, so a widened display doesn't change the score.
    matched, lost, added = multichain.diff_pairs(
        scoped_ref_pairs, model_nested + model_crossing
    )

    def _inf_block(key):
        m = inf_metrics.get(key) or {}
        return {k: m.get(k) for k in ("inf", "ppv", "sty", "tp", "fp", "fn")}

    metrics_payload = {
        "structure_id": structure_id,
        "model_id": model_id,
        "model_simulated": model_is_simulated,
        "model_own_layout": bool(model_own_layout),
        "chains": result.order,
        "inf": {k: _inf_block(k) for k in ("wc", "nwc", "all")},
        "superpose_rmsd": superpose_rmsd,
        "superpose_n_atoms": superpose_n,
        "diff": {
            "matched": len(matched),
            "lost": len(lost),
            "added": len(added),
        },
    }
    (viewer_dir / "metrics.json").write_text(
        json.dumps(metrics_payload, indent=2) + "\n"
    )

    if not quiet:
        rprint(
            f"[cyan]Base-pair diff (reference vs {model_id}): "
            f"{len(matched)} matched, {len(lost)} missing in model, "
            f"{len(added)} model-only[/cyan]"
        )

        def _fmt(metric):
            val = metric.get("inf")
            return "n/a" if val is None else f"{val:.3f}"

        rprint(
            "[cyan]INF (interaction network fidelity): "
            f"WC {_fmt(inf_metrics['wc'])}, "
            f"non-WC {_fmt(inf_metrics['nwc'])}, "
            f"all {_fmt(inf_metrics['all'])}[/cyan]"
        )
        rprint(f"[green]Compare viewer ready: {html_path.resolve()}[/green]")
        rprint(
            "[dim]Serve over HTTP, e.g.:\n"
            f"  python3 -m http.server -d {viewer_dir.resolve()} 8000\n"
            "then open http://localhost:8000/[/dim]"
        )


@cli.command()
@click.argument("pdb-input", type=click.STRING)
@click.argument("output-folder", type=click.Path())
@click.option(
    "--basepairs",
    type=click.Choice(["auto", "rnaview", "fr3d", "cif"]),
    default="auto",
    help=(
        "Tool for base pair extraction (default: auto = prefer FR3D). "
        "'cif' reads pairs from the mmCIF's own DNATCO/NDB annotation, no "
        "FR3D run."
    ),
)
@click.option(
    "--format",
    "structure_format",
    type=click.Choice(["auto", "pdb", "cif"]),
    default="auto",
    help="Preferred structure format to download (default: auto = prefer PDB)",
)
@click.option(
    "--chain",
    type=str,
    default=None,
    help="Specific chain ID to extract (default: first RNA chain)",
)
@click.option(
    "--pseudoknots/--no-pseudoknots",
    default=True,
    help="Include pseudoknots in structure (FR3D only, uses Aa/Bb notation, default: on)",
)
@click.option(
    "--rnapuzzler",
    "rnapuzzler_flag",
    default=False,
    is_flag=True,
    help="Use RNApuzzler for overlap-free layout (ViennaRNA, templatefree only)",
)
@click.option(
    "--mode",
    type=click.Choice(["auto", "templated", "templatefree"]),
    default="auto",
    help=(
        "Layout mode. 'auto' (default) tries the template search first and "
        "falls back to templatefree if no template matches. 'templated' runs "
        "the full R2DT template search (CRW/RiboVision/Rfam/GtRNAdb/RNase P/"
        "tmRNA) and fails if nothing matches. 'templatefree' uses the "
        "FR3D-derived dot-bracket with R2R/RNApuzzler/RNArtist."
    ),
)
@click.option(
    "--all-chains",
    "all_chains",
    is_flag=True,
    default=False,
    help=(
        "Combine every RNA chain into one diagram, auto-ordering the "
        "concatenation to minimise crossing inter-chain pairs (mmCIF only, "
        "templatefree)."
    ),
)
@click.option(
    "--chains",
    "chains",
    type=str,
    default=None,
    help=(
        "Combine the listed RNA chains into one diagram, in the given "
        "concatenation order (e.g. 'A,B'). Implies multi-chain templatefree."
    ),
)
@click.option(
    "--simulate-model",
    "simulate_model",
    is_flag=True,
    default=False,
    help=(
        "TESTING: also emit a <id>.model.svg by randomly perturbing the "
        "reference base pairs, to preview the reference/model diff without a "
        "real predicted structure (multi-chain only)."
    ),
)
@click.option(
    "--simulate-seed",
    "simulate_seed",
    type=int,
    default=2,
    help="Seed for --simulate-model perturbation (default: 2).",
)
@click.option(
    "--compare",
    "compare",
    is_flag=True,
    default=False,
    help=(
        "Emit an interactive 3-panel viewer/ page: reference 2D, model 2D, "
        "and a shared Mol* 3D (multi-chain only). Without --model the model "
        "panel is a simulated perturbation of the reference."
    ),
)
@click.option(
    "--model",
    "model_file",
    type=click.Path(exists=True),
    default=None,
    help=(
        "Predicted model structure (.pdb/.cif) to compare against this "
        "structure as reference. Must share the reference's sequence in the "
        "same chain order. Implies --compare."
    ),
)
@click.option(
    "--model-chains",
    "model_chains",
    type=str,
    default=None,
    help="RNA chains of --model, in order matching --chains (default: all).",
)
@click.option("--quiet", "-q", is_flag=True, default=False)
@click.pass_context
# pylint: disable=too-many-arguments,too-many-branches,too-many-statements,too-many-locals
# pylint: disable=too-many-positional-arguments,too-many-return-statements
def pdb(
    ctx,
    pdb_input,
    output_folder,
    basepairs,
    structure_format,
    chain,
    pseudoknots,
    rnapuzzler_flag,
    mode,
    all_chains,
    chains,
    simulate_model,
    simulate_seed,
    compare,
    model_file,
    model_chains,
    quiet,
):
    """
    Generate R2DT diagram from a PDB structure.

    Accepts either a PDB ID (downloads from RCSB) or a local PDB/mmCIF file.
    Supports gzip-compressed files (.pdb.gz, .cif.gz).

    Examples:

        r2dt.py pdb 1S72 output/

        r2dt.py pdb 9FN3 output/ --basepairs fr3d

        r2dt.py pdb 1EHZ output/ --chain A

        r2dt.py pdb ./structures/my_rna.cif output/

        r2dt.py pdb ./structures/my_rna.pdb.gz output/ --basepairs rnaview
    """
    # Check if input is a local file or PDB ID
    is_local_file = pdb_fetch.is_local_structure_file(pdb_input)

    if is_local_file:
        # Validate local file
        file_path = Path(pdb_input)
        is_valid, detected_format, error_msg = pdb_fetch.validate_structure_file(
            file_path
        )
        if not is_valid:
            rprint(f"[red]Error: {error_msg}[/red]")
            return

        # Derive structure ID from filename
        structure_id = file_path.stem
        if structure_id.endswith((".pdb", ".cif")):
            structure_id = structure_id.rsplit(".", 1)[0]

        actual_format = detected_format

        if not quiet:
            rprint(shared.get_r2dt_version_header())
            rprint(f"[bold]Processing local structure file: {file_path}[/bold]")
            rprint(f"Detected format: {actual_format}")
    else:
        # Treat as PDB ID
        pdb_id = pdb_input
        if not pdb_fetch.validate_pdb_id(pdb_id):
            rprint(
                f"[red]Error: '{pdb_input}' is not a valid PDB ID or file path[/red]"
            )
            rprint("PDB IDs should be 4+ alphanumeric characters (e.g., 1S72)")
            rprint("File paths should end with .pdb, .cif, .pdb.gz, or .cif.gz")
            return

        if not quiet:
            rprint(shared.get_r2dt_version_header())
            rprint(f"[bold]Processing PDB structure: {pdb_id}[/bold]")

        structure_id = pdb_id

    # Create output directory
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)

    if is_local_file:
        # Use local file directly (in-place reference)
        pass  # file_path and actual_format already set above
    else:
        # Download from RCSB
        downloads_dir = output_path / "downloads"
        downloads_dir.mkdir(exist_ok=True)

        # Determine preferred format based on basepairs tool
        if structure_format == "auto":
            # If user wants fr3d/cif, prefer CIF (FR3D works best with CIF; the
            # cif source reads CIF-only annotation). rnaview must use PDB.
            if basepairs in ("fr3d", "cif"):
                prefer_format = "cif"
            else:
                prefer_format = "pdb"
        else:
            prefer_format = structure_format

        if not quiet:
            rprint(f"Downloading structure from RCSB (prefer {prefer_format})...")

        file_path, actual_format = pdb_fetch.download_structure(
            pdb_id, downloads_dir, prefer_format=prefer_format
        )

        if not file_path:
            rprint(f"[red]Error: Could not download structure {pdb_id} from RCSB[/red]")
            return

        if not quiet:
            rprint(f"Downloaded: {file_path} (format: {actual_format})")

    # --- Multi-chain combined diagram (experimental) ---
    if all_chains or chains:
        _run_multichain_pdb(
            ctx,
            file_path=file_path,
            actual_format=actual_format,
            structure_id=structure_id,
            output_path=output_path,
            chains=chains,
            rnapuzzler_flag=rnapuzzler_flag,
            simulate_model=simulate_model or compare or bool(model_file),
            simulate_seed=simulate_seed,
            compare=compare or bool(model_file),
            model_file=model_file,
            model_chains=model_chains,
            quiet=quiet,
        )
        return

    # Determine which basepairs tool to use
    if basepairs == "auto":
        # Prefer FR3D: it supports pseudoknots and handles both PDB and CIF
        use_basepairs = "fr3d"
    elif basepairs == "rnaview" and actual_format == "cif":
        rprint(
            "[red]Error: RNAView cannot process mmCIF files. "
            "Use --basepairs fr3d or --format pdb[/red]"
        )
        return
    elif basepairs == "cif":
        if actual_format != "cif":
            rprint(
                "[red]Error: --basepairs cif needs an mmCIF input "
                "(use --format cif or provide a .cif file)[/red]"
            )
            return
        if not cif_basepairs.has_annotation(file_path):
            rprint(
                "[red]Error: this mmCIF has no base-pair annotation "
                "(_ndb_base_pair_list / _ndb_base_pair_annotation). Use a "
                "DNATCO/NDB-annotated CIF, or --basepairs fr3d[/red]"
            )
            return
        use_basepairs = "cif"
    else:
        use_basepairs = basepairs

    if not quiet:
        rprint(f"Using {use_basepairs} for base pair extraction...")

    # Pre-detect RNA chain so FR3D and RNAview use the same one
    if not chain and actual_format == "pdb":
        _, _, detected_chain = fr3d_utils.get_full_sequence(str(file_path), None)
        if detected_chain:
            chain = detected_chain
            if not quiet:
                rprint(f"Auto-detected RNA chain: {chain}")

    # Extract secondary structure
    extraction_dir = output_path / "extraction"
    extraction_dir.mkdir(exist_ok=True)

    if use_basepairs == "rnaview":
        # Use existing rnaview module
        sequence, dot_bracket = _extract_with_rnaview(str(file_path), chain, quiet)
    elif use_basepairs == "cif":
        # Read pairs from the CIF's own annotation; no FR3D run.
        sequence, dot_bracket = cif_basepairs.get_secondary_structure_cif(
            str(file_path),
            str(extraction_dir),
            structure_id=structure_id,
            chain_id=chain,
            include_pseudoknots=pseudoknots,
            quiet=quiet,
        )
    else:
        # Use FR3D
        sequence, dot_bracket = fr3d_utils.get_secondary_structure_fr3d(
            str(file_path),
            str(extraction_dir),
            chain_id=chain,
            include_pseudoknots=pseudoknots,
            quiet=quiet,
        )

        # If FR3D returned a sequence but no base pairs (crash/failure),
        # fall back to RNAview for PDB files
        if (
            sequence
            and dot_bracket
            and dot_bracket.count("(") == 0
            and actual_format == "pdb"
        ):
            if not quiet:
                rprint(
                    "[yellow]FR3D found no base pairs, "
                    "trying RNAView as fallback...[/yellow]"
                )
            rv_seq, rv_db = _extract_with_rnaview(str(file_path), chain, quiet)
            if rv_seq and rv_db and rv_db.count("(") > 0:
                sequence, dot_bracket = rv_seq, rv_db
                if not quiet:
                    rprint(
                        f"[green]RNAView recovered {rv_db.count('(')} "
                        f"base pairs[/green]"
                    )

    if not sequence or not dot_bracket:
        rprint("[red]Error: Could not extract secondary structure[/red]")
        return

    if not quiet:
        rprint(f"Sequence length: {len(sequence)}")
        rprint(f"Base pairs: {dot_bracket.count('(')}")

    # --- Full-sequence expansion (missing residues) ---
    resolved_mask = None
    full_sequence, mask, _used_chain = fr3d_utils.get_full_sequence(
        str(file_path), chain
    )
    if full_sequence and len(full_sequence) > len(sequence):
        n_missing = len(full_sequence) - len(sequence)
        if not quiet:
            rprint(
                f"Full deposited sequence: {len(full_sequence)} nt "
                f"({n_missing} unresolved)"
            )
        try:
            dot_bracket = fr3d_utils.remap_dot_bracket(dot_bracket, mask)
            sequence = full_sequence
            resolved_mask = mask
        except ValueError as exc:
            if not quiet:
                rprint(
                    f"[yellow]Could not map to full sequence: {exc}. "
                    f"Using resolved sequence only.[/yellow]"
                )
    elif not quiet and full_sequence:
        rprint("No missing residues detected")

    # Write FASTA file for R2DT (3-line: seq + dot-bracket, for templatefree)
    fasta_path = output_path / f"{structure_id}.fasta"
    with open(fasta_path, "w") as f:
        f.write(f">{structure_id}\n")
        f.write(f"{sequence}\n")
        f.write(f"{dot_bracket}\n")

    if not quiet:
        rprint(f"Created FASTA: {fasta_path}")

    results_folder = output_path / "results"

    def run_templatefree():
        # Hand the FR3D dot-bracket to R2R / RNApuzzler / RNArtist.
        if not quiet:
            rprint("Generating 2D diagram with R2DT (templatefree mode)...")
        ctx.invoke(
            templatefree,
            fasta_input=str(fasta_path),
            output_folder=str(results_folder),
            rnartist=False,
            rscape=not rnapuzzler_flag,
            rnapuzzler_flag=rnapuzzler_flag,
            quiet=quiet,
        )

    def run_templated():
        # Feed a sequence-only fasta into `draw`, which picks a template
        # (CRW / RiboVision / Rfam / GtRNAdb / RNase P / tmRNA / Rfam-tRNA)
        # and lays the diagram out via Traveler. Returns the matched
        # template id, or None if nothing matched.
        if not quiet:
            rprint("Generating 2D diagram with R2DT (templated mode)...")
        draw_fasta = output_path / f"{structure_id}.draw.fasta"
        with open(draw_fasta, "w") as f:
            f.write(f">{structure_id}\n")
            f.write(f"{sequence}\n")
        ctx.invoke(
            draw,
            fasta_input=str(draw_fasta),
            output_folder=str(results_folder),
            quiet=quiet,
        )
        # `draw` names its outputs ``<structure_id>-<template_id>.colored.*``.
        # Downstream code (grey-out, viewer-export) keys on the plain
        # ``<structure_id>`` basename, so collapse the template suffix here.
        return _rename_templated_outputs(results_folder, structure_id)

    if mode == "templatefree":
        run_templatefree()
    else:
        # 'templated' and 'auto' both try the template search first.
        matched_template = run_templated()
        if matched_template is not None:
            if not quiet:
                rprint(f"[green]Matched template: {matched_template}[/green]")
        elif mode == "templated":
            rprint("[red]No template matched this structure in templated mode.[/red]")
            rprint(
                "[yellow]Try --mode templatefree, or --mode auto, which "
                "falls back to the FR3D-derived layout automatically.[/yellow]"
            )
            return
        else:
            # auto: no template matched -> fall back to templatefree. Clear
            # the partial `draw` output first so it can't shadow the
            # templatefree results.
            if not quiet:
                rprint(
                    "[yellow]No template matched; falling back to "
                    "templatefree layout.[/yellow]"
                )
            if results_folder.exists():
                shutil.rmtree(results_folder)
            run_templatefree()

    # --- Post-process: grey out unresolved nucleotides ---
    if resolved_mask is not None:
        n_greyed = pdb_post.grey_out_svg_files(
            str(results_folder), structure_id, resolved_mask, quiet=quiet
        )
        if not quiet and n_greyed:
            rprint(
                f"[dim]Greyed out unresolved nucleotides in "
                f"{n_greyed} SVG file(s)[/dim]"
            )

    # Report success
    svg_path = results_folder / "results" / "svg" / f"{structure_id}.svg"
    if svg_path.exists():
        if not quiet:
            rprint(f"[green]Success! SVG created: {svg_path}[/green]")
    else:
        if not quiet:
            rprint(
                "[yellow]Diagram generation completed. Check output folder.[/yellow]"
            )


@cli.command("pdb_2d_3d")
@click.argument("pdb-input", type=click.STRING)
@click.argument("output-folder", type=click.Path())
@click.option(
    "--basepairs",
    type=click.Choice(["auto", "rnaview", "fr3d", "cif"]),
    default="auto",
)
@click.option(
    "--format",
    "structure_format",
    type=click.Choice(["auto", "pdb", "cif"]),
    default="auto",
)
@click.option("--chain", type=str, default=None)
@click.option("--pseudoknots/--no-pseudoknots", default=True)
@click.option("--rnapuzzler", "rnapuzzler_flag", default=False, is_flag=True)
@click.option(
    "--mode",
    type=click.Choice(["auto", "templated", "templatefree"]),
    default="auto",
)
@click.option("--quiet", "-q", is_flag=True, default=False)
@click.pass_context
# pylint: disable=too-many-arguments,too-many-locals,too-many-positional-arguments
# pylint: disable=too-many-statements,too-many-branches
def pdb_2d_3d(
    ctx,
    pdb_input,
    output_folder,
    basepairs,
    structure_format,
    chain,
    pseudoknots,
    rnapuzzler_flag,
    mode,
    quiet,
):
    """
    Generate an interactive 2D+3D viewer page from a PDB structure.

    Runs the same pipeline as ``pdb`` and additionally writes a
    ``viewer/`` folder containing ``index.html``, the two JSON blobs
    consumed by pdb-rna-viewer, the structure file, and the viewer
    assets.  Opening ``index.html`` in a browser shows the linked 2D
    diagram + 3D molstar view.

    Examples:

        r2dt.py pdb_2d_3d 1Y26 output/

        r2dt.py pdb_2d_3d ./my_rna.cif output/ --basepairs fr3d
    """
    # 1) Run the existing pdb pipeline.
    ctx.invoke(
        pdb,
        pdb_input=pdb_input,
        output_folder=output_folder,
        basepairs=basepairs,
        structure_format=structure_format,
        chain=chain,
        pseudoknots=pseudoknots,
        rnapuzzler_flag=rnapuzzler_flag,
        mode=mode,
        quiet=quiet,
    )

    # 2) Resolve the file paths the pdb command wrote.
    output_path = Path(output_folder)
    is_local_file = pdb_fetch.is_local_structure_file(pdb_input)
    if is_local_file:
        file_path = Path(pdb_input)
        structure_id = file_path.stem
        if structure_id.endswith((".pdb", ".cif")):
            structure_id = structure_id.rsplit(".", 1)[0]
        _, actual_format, _ = pdb_fetch.validate_structure_file(file_path)
        source_structure_path = file_path
    else:
        structure_id = pdb_input
        downloads = output_path / "downloads"
        # The pdb command preserves the original-case PDB ID in the
        # downloaded filename; lowercasing it would break the unit_id
        # key match against the FR3D basepair file.
        cands = (
            list(downloads.glob(f"{structure_id}.cif"))
            + list(downloads.glob(f"{structure_id}.pdb"))
            + list(downloads.glob(f"{structure_id}.*"))
        )
        cands = [c for c in cands if c.is_file()]
        if not cands:
            rprint("[red]Viewer step: cannot locate downloaded structure[/red]")
            return
        source_structure_path = cands[0]
        actual_format = source_structure_path.suffix.lstrip(".")

    colored_json = (
        output_path / "results" / "results" / "json" / f"{structure_id}.colored.json"
    )
    colored_svg = (
        output_path / "results" / "results" / "svg" / f"{structure_id}.colored.svg"
    )
    basepair_txt = output_path / "extraction" / f"{structure_id}_basepair.txt"
    # If `pdb` bailed (e.g. templated mode with no template match) there
    # is no colored SVG/JSON to build a viewer from -- abort quietly; the
    # pdb step has already explained why.
    if not colored_json.exists():
        return

    # 3) Re-derive resolved_mask and unit_id_to_position so we can write
    # both the apiData mapping and the FR3D label remap.  These calls
    # are the same the pdb command makes internally.
    if basepairs == "cif":
        # CIF source: read sequence/positions from the vendored script (no
        # fr3d-python), keyed identically to the basepair file written above.
        cif_chain = cif_basepairs.resolve_chain(str(source_structure_path), chain)
        _, resolved_mask, _ = fr3d_utils.get_full_sequence(
            str(source_structure_path), cif_chain
        )
        _, unit_id_to_position, used_chain = cif_basepairs.read_sequence_and_positions(
            str(source_structure_path), cif_chain, quiet=quiet
        )
    else:
        _, resolved_mask, used_chain = fr3d_utils.get_full_sequence(
            str(source_structure_path), chain
        )
        if str(source_structure_path).lower().endswith(".cif"):
            _, unit_id_to_position = fr3d_utils.extract_sequence_from_cif(
                str(source_structure_path), used_chain or chain, quiet=quiet
            )
        else:
            _, unit_id_to_position = fr3d_utils.extract_sequence_from_pdb(
                str(source_structure_path), used_chain or chain, quiet=quiet
            )
    if not resolved_mask:
        resolved_mask = None

    # Fall back to deriving the chain from a FR3D unit_id when neither the
    # user nor the upstream auto-detect supplied one. Without this the
    # viewer ends up with auth_asym_id="" and molstar's visual.select()
    # matches no residue.
    effective_chain = used_chain or chain
    if not effective_chain and unit_id_to_position:
        for unit in unit_id_to_position:
            parts = unit.split("|")
            if len(parts) >= 3 and parts[2]:
                effective_chain = parts[2]
                break

    # 4) Build the JSON blobs.
    colored = json.loads(colored_json.read_text())
    n_full = sum(
        1
        for nuc in colored["rnaComplexes"][0]["rnaMolecules"][0]["sequence"]
        if nuc.get("residueName") not in ("5'", "3'")
        and len(nuc.get("residueName", "")) == 1
    )
    api_data = viewer_export.build_api_data(
        colored_json,
        structure_id=structure_id,
        chain_id=effective_chain,
        resolved_mask=resolved_mask,
        unit_id_to_position=unit_id_to_position,
        colored_svg_path=colored_svg if colored_svg.exists() else None,
    )
    fr3d_data = viewer_export.build_fr3d_data(
        basepair_txt,
        structure_id=structure_id,
        chain_id=effective_chain,
        unit_id_to_position=unit_id_to_position or {},
        resolved_mask=resolved_mask,
        n_full=n_full,
    )

    # 5) Lay out the viewer folder.
    viewer_dir = output_path / "viewer"
    viewer_dir.mkdir(exist_ok=True)
    (viewer_dir / "api.json").write_text(json.dumps(api_data))
    (viewer_dir / "fr3d.json").write_text(json.dumps(fr3d_data))
    lbn_data = lbn_export.build_lbn_data(api_data, fr3d_data)
    (viewer_dir / "lbn.json").write_text(json.dumps(lbn_data))
    structure_dest_name = f"{structure_id}.{actual_format}"
    shutil.copyfile(source_structure_path, viewer_dir / structure_dest_name)

    # Copy the vendored pdb-rna-viewer build files next to index.html.
    _copy_viewer_assets(viewer_dir)

    html_path = viewer_html.render(
        viewer_dir,
        structure_id=structure_id,
        chain_id=effective_chain,
        structure_filename=structure_dest_name,
        structure_format=actual_format,
        annotation_source=viewer_html.ANNOTATION_SOURCE_HTML.get(
            basepairs if basepairs != "auto" else "fr3d"
        ),
    )

    if not quiet:
        rprint(f"[dim]Viewer chain: {effective_chain or '(none)'}[/dim]")
        # The viewer fetches api.json / fr3d.json / the structure file via
        # relative URLs, which browsers block over file://. Tell the user
        # to serve the folder instead of double-clicking the HTML.
        rprint(f"[green]Viewer ready: {html_path.resolve()}[/green]")
        rprint(
            "[dim]Serve it over HTTP, e.g.:\n"
            f"  python3 -m http.server -d {viewer_dir.resolve()} 8000\n"
            "then open http://localhost:8000/[/dim]"
        )


def _copy_viewer_assets(viewer_dir: Path) -> None:
    """Copy the vendored viewer assets next to ``index.html``.

    The pdb-rna-viewer compiled bundle isn't on a CDN, isn't on npm, and
    the GitHub release downloads are served with ``application/octet-stream``
    which browsers refuse to load as a stylesheet. So we vendor it (plus
    the ``r2dt-2d-3d-viewer.js`` interaction glue) under ``data/viewer/`` in the R2DT
    repo (Apache-2.0) and copy it into each output folder.
    """
    src = Path(__file__).resolve().parent / "data" / "viewer"
    wanted = (
        viewer_html.VIEWER_PLUGIN_FILENAME,
        viewer_html.VIEWER_CSS_FILENAME,
        viewer_html.R2DT_CSS_FILENAME,
        viewer_html.VIEWER_JS_FILENAME,
    )
    missing = [name for name in wanted if not (src / name).is_file()]
    if missing:
        raise click.ClickException(
            f"Missing vendored viewer assets in {src}: {', '.join(missing)}. "
            "The R2DT checkout looks incomplete."
        )
    for name in wanted:
        shutil.copyfile(src / name, viewer_dir / name)


def _rename_templated_outputs(results_folder: Path, structure_id: str):
    """Collapse ``<structure_id>-<template_id>.<ext>`` filenames produced
    by ``draw`` into plain ``<structure_id>.<ext>`` so the rest of the
    ``pdb`` / ``pdb_2d_3d`` pipeline (grey-out, viewer-export) doesn't
    need to know which template won.

    Returns the matched template id (the trailing portion after the
    structure id) if a templated SVG was found, otherwise ``None``.
    """
    svg_dir = results_folder / "results" / "svg"
    candidates = list(svg_dir.glob(f"{structure_id}-*.colored.svg"))
    if not candidates:
        return None
    # Take the first match -- there should be only one per structure.
    colored_svg = candidates[0]
    full_stem = colored_svg.name[: -len(".colored.svg")]
    template_id = full_stem[len(structure_id) + 1 :]

    suffix_dirs = {
        "results/svg": [".colored.svg", ".enriched.svg"],
        "results/thumbnail": [".thumbnail.svg"],
        "results/json": [".colored.json"],
        "results/fasta": [".fasta"],
    }
    for subdir, suffixes in suffix_dirs.items():
        d = results_folder / subdir
        if not d.is_dir():
            continue
        for suffix in suffixes:
            src = d / f"{full_stem}{suffix}"
            if src.exists():
                src.replace(d / f"{structure_id}{suffix}")
    return template_id


def _extract_with_rnaview(pdb_file: str, chain_id=None, quiet=False):
    """
    Extract secondary structure using RNAView.

    Supports gzip-compressed PDB files (.pdb.gz).

    Args:
        pdb_file: Path to PDB file (may be gzip-compressed).
        chain_id: Optional chain ID. If None, uses first chain only.
        quiet: If True, suppress verbose output.

    Returns:
        Tuple of (sequence, dot_bracket) or (None, None).
    """
    try:
        # Use DecompressedStructureFile to handle .gz files
        # RNAView requires uncompressed file on disk
        with pdb_fetch.DecompressedStructureFile(Path(pdb_file)) as decompressed_path:
            # Extract sequence using rnaview module
            # If no chain specified, use first chain only (consistent with FR3D behavior)
            sequence = rnaview_utils.extract_sequence(
                str(decompressed_path), chain_id=chain_id, quiet=quiet
            )

            if not sequence:
                return None, None

            # Run RNAView
            rnaview_output = rnaview_utils.run_rnaview(str(decompressed_path))

            # Parse output to dot-bracket
            dot_bracket = rnaview_utils.parse_rnaview_output(
                rnaview_output, sequence, quiet
            )

            # Clean up temporary files created by RNAView
            rnaview_utils.cleanup_rnaview_files(str(decompressed_path))

            return sequence, dot_bracket

    except Exception as e:  # pylint: disable=broad-exception-caught
        if not quiet:
            print(f"Error in RNAView extraction: {e}")
        return None, None


@cli.command("workstation")
@click.option(
    "--workspace",
    type=click.Path(),
    default=None,
    help="Local cache directory (default: ~/.r2dt-workstation).",
)
@click.option("--port", default=8765, show_default=True, type=int)
@click.option(
    "--bind",
    default="127.0.0.1",
    show_default=True,
    help="Bind address. Use 0.0.0.0 only inside Docker with -p 127.0.0.1:PORT:PORT.",
)
@click.option(
    "--docker-image",
    default="rnacentral/r2dt:latest",
    show_default=True,
    help="Image used when the server runs on the host and spawns job containers.",
)
def workstation(workspace, port, bind, docker_image):
    """
    Start the private local R2DT workstation (web UI).

    Homepage plus mode dashboards for 2D, 2D+3D, compare, and alignments.
    Requires Docker. Prefer ``just workstation``, which publishes the port
    to 127.0.0.1 only and mounts ~/.r2dt-workstation.
    """
    ws = (
        Path(workspace).expanduser() if workspace else Path.home() / ".r2dt-workstation"
    )
    repo_root = Path(__file__).resolve().parent
    workstation_mod.run_server(
        workspace=ws,
        repo_root=repo_root,
        host=bind,
        port=port,
        docker_image=docker_image,
    )


if __name__ == "__main__":
    cli()
