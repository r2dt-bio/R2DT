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
import tarfile
import time
import unittest
from pathlib import Path

import click  # pylint: disable=import-error
from rich import print as rprint

from tests import tests
from utils import config, core
from utils import generate_cm_library as gcl
from utils import generate_model_info as gmi
from utils import gtrnadb
from utils import list_models as lm
from utils import r2r, rfam
from utils import rna2djsonschema as r2djs
from utils import shared
from utils.runner import runner


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
    if not os.path.exists(config.CM_LIBRARY):
        os.makedirs(config.CM_LIBRARY)
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
    destination_dir = os.path.join(config.CM_LIBRARY, "crw")

    if os.path.exists(source_dir):
        shutil.move(source_dir, destination_dir)

    rprint("Generating CRW modelinfo file")
    gmi.generate_model_info(cm_library=config.CRW_CM_LIBRARY)


@cli.command()
def setup_rfam():
    """
    Re-generate Rfam templates from scratch.
    """
    rprint(shared.get_r2dt_version_header())
    rfam.setup()
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
    runner.run(f"esl-sfetch --index {output_filename}")


def is_templatefree(fasta_input):
    """Check if the input file is a valid fasta file
    with an additional line specifying secondary structure
    in dot bracket format (pseudoknots allowed)."""
    with open(fasta_input) as f_in:
        lines = [line for line in f_in.readlines() if line.strip()]
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
# pylint: disable-next=too-many-arguments, too-many-locals, too-many-statements, too-many-branches
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
    if not quiet:
        rprint(shared.get_r2dt_version_header())

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

    hits = set()
    subset_fasta = os.path.join(output_folder, "subset.fasta")
    shutil.copy(fasta_input, subset_fasta)
    runner.run(f"esl-sfetch --index {subset_fasta}")

    # Rfam
    if not quiet:
        rprint(f"Analysing {len(all_seq_ids)} sequences with [yellow]Rfam[/yellow]")
    with Timer("Rfam", quiet):
        with open(
            shared.get_ribotyper_output(
                fasta_input,
                rfam_output,
                os.path.join(config.CM_LIBRARY, "rfam"),
                skip_ribovore_filters,
            ),
        ) as f_ribotyper:
            for line in f_ribotyper.readlines():
                rnacentral_id, model_id, _ = line.split("\t")
                core.visualise(
                    "rfam",
                    fasta_input,
                    rfam_output,
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

    # RiboVision SSU
    hits = hits.union(get_hits(rfam_output))
    subset = all_seq_ids.difference(hits)
    if subset:
        get_subset_fasta(fasta_input, subset_fasta, subset)
        with Timer("RV SSU", quiet):
            if not quiet:
                rprint(f"Analysing {len(subset)} sequences with RiboVision SSU")
            ctx.invoke(
                ribovision_draw_ssu,
                fasta_input=subset_fasta,
                output_folder=ribovision_ssu_output,
                constraint=constraint,
                exclusion=exclusion,
                fold_type=fold_type,
                quiet=True,
                skip_ribovore_filters=skip_ribovore_filters,
            )

    # CRW
    hits = hits.union(get_hits(ribovision_ssu_output))
    subset = all_seq_ids.difference(hits)
    if subset:
        get_subset_fasta(fasta_input, subset_fasta, subset)
        with Timer("CRW", quiet):
            if not quiet:
                rprint(f"Analysing {len(subset)} sequences with CRW")
            ctx.invoke(
                rrna_draw,
                fasta_input=subset_fasta,
                output_folder=crw_output,
                constraint=constraint,
                exclusion=exclusion,
                fold_type=fold_type,
                quiet=True,
                skip_ribovore_filters=skip_ribovore_filters,
            )

    # RiboVision LSU
    hits = hits.union(get_hits(crw_output))
    subset = all_seq_ids.difference(hits)
    if subset:
        get_subset_fasta(fasta_input, subset_fasta, subset)
        with Timer("LSU", quiet):
            if not quiet:
                rprint(f"Analysing {len(subset)} sequences with RiboVision LSU")
            ctx.invoke(
                ribovision_draw_lsu,
                fasta_input=subset_fasta,
                output_folder=ribovision_lsu_output,
                constraint=constraint,
                exclusion=exclusion,
                fold_type=fold_type,
                quiet=True,
                skip_ribovore_filters=skip_ribovore_filters,
            )

    # RNAse P
    hits = hits.union(get_hits(ribovision_lsu_output))
    subset = all_seq_ids.difference(hits)
    if subset:
        get_subset_fasta(fasta_input, subset_fasta, subset)
        with Timer("RNAse P", quiet):
            if not quiet:
                rprint(f"Analysing {len(subset)} sequences with RNAse P models")
            ctx.invoke(
                rnasep_draw,
                fasta_input=subset_fasta,
                output_folder=rnasep_output,
                constraint=constraint,
                exclusion=exclusion,
                fold_type=fold_type,
                quiet=True,
                skip_ribovore_filters=skip_ribovore_filters,
            )

    # GtRNAdb
    hits = hits.union(get_hits(rnasep_output))
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
    ]
    for folder in result_folders:
        organise_results(folder, output_folder)
    organise_metadata(output_folder, result_folders)

    # clean up
    os.system(f"rm {output_folder}/subset*")
    os.system(f"rm -f {fasta_input}.ssi")


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
    # pylint: disable=too-many-arguments
    if not quiet:
        rprint(shared.get_r2dt_version_header())
    os.makedirs(output_folder, exist_ok=True)

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
    # pylint: disable=too-many-arguments
    if not quiet:
        rprint(shared.get_r2dt_version_header())
    os.makedirs(output_folder, exist_ok=True)
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
    # pylint: disable=too-many-arguments
    if not quiet:
        rprint(shared.get_r2dt_version_header())
    os.makedirs(output_folder, exist_ok=True)
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
    # pylint: disable=too-many-arguments
    if not quiet:
        rprint(shared.get_r2dt_version_header())
    os.makedirs(output_folder, exist_ok=True)
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
    # pylint: disable=too-many-arguments
    if not quiet:
        rprint(shared.get_r2dt_version_header())
    os.makedirs(output_folder, exist_ok=True)
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
):
    """
    Visualise sequences using the Rfam/R-scape consensus structure as template.

    RFAM_ACCESSION - Rfam family to process (RF00001, RF00002 etc)
    """
    # pylint: disable=too-many-arguments
    if not quiet:
        rprint(shared.get_r2dt_version_header())
        rprint(rfam_acc)
    if rnartist:
        template_type = "rnartist"
    else:
        template_type = "rscape"
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
    # pylint: disable=too-many-arguments, too-many-locals
    if not quiet:
        rprint(shared.get_r2dt_version_header())
    model_type = lm.get_model_type(model_id)
    if not model_type:
        rprint("Error: Model not found. Please check model_id")
        return
    if not quiet:
        rprint(
            f"Visualising sequence {seq_id} using the {model_id} model from {model_type}"
        )
    runner.run(f"esl-sfetch --index {fasta_input}")

    output = os.path.join(output_folder, model_type.replace("_", "-"))

    if model_type in ["crw", "rfam", "rnasep", "local_data"]:
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
        "local_data": "Local data",
    }
    with open(
        os.path.join(metadata_folder, "metadata.tsv"), "a", encoding="utf-8"
    ) as f_out:
        line = f"{seq_id}\t{model_id}\t{label_mapping[model_type]}\n"
        f_out.write(line)


@cli.command()
@click.argument("fasta-input", type=click.Path())
@click.argument("output-folder", type=click.Path())
@click.option("--quiet", "-q", is_flag=True, default=False)
def templatefree(fasta_input, output_folder, quiet):
    """
    Run template-free visualisation using R2R to generate a layout.
    """
    if not quiet:
        rprint(shared.get_r2dt_version_header())
    results_folder = os.path.join(output_folder, "results")
    r2r_folder = os.path.join(output_folder, "r2r")
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(results_folder, exist_ok=True)
    os.makedirs(r2r_folder, exist_ok=True)
    seq_id, sequence, structure = r2r.parse_fasta(fasta_input)
    r2r.generate_r2r_input_file(sequence, structure, r2r_folder)
    r2r_svg = r2r.run_r2r(r2r_folder)
    rscape_one_line_svg = rfam.convert_rscape_svg_to_one_line(r2r_svg, r2r_folder)
    rfam.convert_rscape_svg_to_traveler(rscape_one_line_svg, r2r_folder)
    r2r.run_traveler(fasta_input, r2r_folder, seq_id)
    organise_results(r2r_folder, output_folder)
    tsv_folder = os.path.join(results_folder, "tsv")
    os.makedirs(tsv_folder, exist_ok=True)
    with open(os.path.join(tsv_folder, "metadata.tsv"), "w") as f_out:
        f_out.write(f"{seq_id}\tR2R\tR2R\n")
    shutil.copyfile(
        fasta_input,
        os.path.join(results_folder, "fasta", f"{seq_id}.fasta"),
    )


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


if __name__ == "__main__":
    cli()
