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
from utils import config, core
from utils import fr3d as fr3d_utils
from utils import generate_cm_library as gcl
from utils import generate_model_info as gmi
from utils import gtrnadb
from utils import list_models as lm
from utils import pdb_fetch, r2r, rfam
from utils import rna2djsonschema as r2djs
from utils import rnaview as rnaview_utils
from utils import shared
from utils import stockholm as stockholm_utils
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
    Run both R2R and RNArtist, return the layout with fewer overlaps.
    Prefers R2R if overlap counts are equal or very similar.
    """
    # pylint: disable=import-outside-toplevel,too-many-branches,too-many-statements,too-many-locals
    import tempfile

    fasta_input = shared.sanitise_fasta(fasta_input)
    seq_id, sequence, structure = r2r.parse_fasta(fasta_input)

    output_folder = Path(output_folder)
    results_folder = output_folder / "results"

    # Create temp directories for both layouts
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        r2r_output = tmpdir / "r2r"
        rnartist_output = tmpdir / "rnartist"

        # Run R2R
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

        # Check R2R overlaps first - skip RNArtist if R2R is perfect
        r2r_overlaps = _count_overlaps(r2r_json)

        if r2r_overlaps == 0:
            # R2R has no overlaps, no need to run RNArtist
            if not quiet:
                rprint("[green]R2R overlaps: 0 -> Using R2R[/green]")
            winner = "R2R"
            winner_output = r2r_output
        else:
            # Run RNArtist to compare
            rnartist_folder_work = rnartist_output / "rnartist"
            for folder in [rnartist_output, rnartist_folder_work]:
                folder.mkdir(exist_ok=True, parents=True)

            rnartist_obj = RnaArtist(destination=rnartist_folder_work)
            rnartist_obj.fasta_file = fasta_input
            rnartist_obj.seq_label = seq_id
            rnartist_obj.run(rerun=True)

            # Scale RNArtist coordinates to avoid Traveler numerical issues
            # RNArtist can place nucleotides very close together in tight loops
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

            # Fix font size in RNArtist SVGs (Traveler uses tiny fonts for tight layouts)
            rnartist_svg = rnartist_output / "results" / "svg" / f"{seq_id}.colored.svg"
            if rnartist_svg.exists():
                _fix_svg_font_size(rnartist_svg)

            rnartist_json = (
                rnartist_output / "results" / "json" / f"{seq_id}.colored.json"
            )

            # Count RNArtist overlaps
            rnartist_overlaps = _count_overlaps(rnartist_json)

            # Choose winner (prefer R2R if equal or R2R has fewer)
            if r2r_overlaps <= rnartist_overlaps:
                winner = "R2R"
                winner_output = r2r_output
            else:
                winner = "RNArtist"
                winner_output = rnartist_output

            if not quiet:
                rprint(
                    f"[green]R2R overlaps: {r2r_overlaps}, "
                    f"RNArtist overlaps: {rnartist_overlaps} -> Using {winner}[/green]"
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


@cli.command()
@click.argument("fasta-input", type=click.Path())
@click.argument("output-folder", type=click.Path())
@click.option("--quiet", "-q", is_flag=True, default=False)
@click.option("--rnartist", default=False, is_flag=True)
@click.option("--rscape", default=False, is_flag=True)
@click.option(
    "--auto",
    default=False,
    is_flag=True,
    help="Run both R2R and RNArtist, pick best layout",
)
# pylint: disable=too-many-arguments,too-many-positional-arguments
# pylint: disable=too-many-statements,inconsistent-return-statements
def templatefree(fasta_input, output_folder, rnartist, rscape, auto, quiet):
    """
    Run template-free visualisation using R2R to generate a layout.
    """
    if not quiet:
        rprint(shared.get_r2dt_version_header())

    # Handle auto mode
    if auto:
        return _templatefree_auto(fasta_input, output_folder, quiet)

    if not rnartist and not rscape:
        rscape = True
    if rnartist and rscape:
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
        r2r.run_traveler(fasta_input, r2r_folder, seq_id)
        organise_results(r2r_folder, output_folder)
        tsv_folder = results_folder / "tsv"
        tsv_folder.mkdir(exist_ok=True, parents=True)
        with open(tsv_folder / "metadata.tsv", "w") as f_out:
            f_out.write(f"{seq_id}\tR2R\tR2R\n")
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
    "--auto-repair/--no-auto-repair",
    default=False,
    help="Attempt to auto-repair unbalanced bracket structures",
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
    auto_repair,
):
    """
    Process a Stockholm alignment with named secondary structure regions.

    Parses the Stockholm file, extracts regions marked in #=GC knownSS_names,
    computes RF consensus sequences, and generates template-free visualizations
    for each region. Optionally stitches all outputs into a single SVG.

    Example:
        r2dt.py stockholm alignment.stk output_folder
    """
    # pylint: disable=import-outside-toplevel,redefined-outer-name
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
            # Generate captions from region names
            captions = [r["name"] for r in processed_regions]

            # Stitch
            combined = stitch_module.stitch_svgs(
                panels=panels,
                gap=100,
                glyph_type="break",
                captions=captions,
                monochrome=monochrome,
                show_gap_labels=True,
                outline=True,
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


@cli.command()
@click.argument("pdb-input", type=click.STRING)
@click.argument("output-folder", type=click.Path())
@click.option(
    "--basepairs",
    type=click.Choice(["auto", "rnaview", "fr3d"]),
    default="auto",
    help="Tool for base pair extraction (default: auto)",
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
    "--pseudoknots",
    is_flag=True,
    default=False,
    help="Include pseudoknots in structure (FR3D only, uses Aa/Bb notation)",
)
@click.option("--quiet", "-q", is_flag=True, default=False)
@click.pass_context
# pylint: disable=too-many-arguments,too-many-branches,too-many-statements,too-many-locals
# pylint: disable=too-many-positional-arguments
def pdb(
    ctx,
    pdb_input,
    output_folder,
    basepairs,
    structure_format,
    chain,
    pseudoknots,
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
            # If user wants fr3d, prefer CIF (FR3D works best with CIF)
            # If user wants rnaview, must use PDB
            if basepairs == "fr3d":
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

    # Determine which basepairs tool to use
    if basepairs == "auto":
        if actual_format == "pdb":
            use_basepairs = "rnaview"
        else:
            use_basepairs = "fr3d"
    elif basepairs == "rnaview" and actual_format == "cif":
        rprint(
            "[red]Error: RNAView cannot process mmCIF files. "
            "Use --basepairs fr3d or --format pdb[/red]"
        )
        return
    else:
        use_basepairs = basepairs

    if not quiet:
        rprint(f"Using {use_basepairs} for base pair extraction...")

    # Extract secondary structure
    extraction_dir = output_path / "extraction"
    extraction_dir.mkdir(exist_ok=True)

    if use_basepairs == "rnaview":
        # Use existing rnaview module
        sequence, dot_bracket = _extract_with_rnaview(str(file_path), chain, quiet)
    else:
        # Use FR3D
        sequence, dot_bracket = fr3d_utils.get_secondary_structure_fr3d(
            str(file_path),
            str(extraction_dir),
            chain_id=chain,
            include_pseudoknots=pseudoknots,
            quiet=quiet,
        )

    if not sequence or not dot_bracket:
        rprint("[red]Error: Could not extract secondary structure[/red]")
        return

    if not quiet:
        rprint(f"Sequence length: {len(sequence)}")
        rprint(f"Base pairs: {dot_bracket.count('(')}")

    # Write FASTA file for R2DT
    fasta_path = output_path / f"{structure_id}.fasta"
    with open(fasta_path, "w") as f:
        f.write(f">{structure_id}\n")
        f.write(f"{sequence}\n")
        f.write(f"{dot_bracket}\n")

    if not quiet:
        rprint(f"Created FASTA: {fasta_path}")

    # Run R2DT templatefree
    if not quiet:
        rprint("Generating 2D diagram with R2DT...")

    results_folder = output_path / "results"
    ctx.invoke(
        templatefree,
        fasta_input=str(fasta_path),
        output_folder=str(results_folder),
        rnartist=False,
        rscape=True,
        quiet=quiet,
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


if __name__ == "__main__":
    cli()
