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

import glob
import os

import click

from utils import crw, rfam, ribovision, gtrnadb, config
from utils.generate_model_info import generate_model_info


def get_ribotyper_output(fasta_input, output_folder, cm_library):
    """
    Run ribotyper on the fasta sequences to select the best matching covariance
    model.
    """
    ribotyper_long_out = os.path.join(output_folder, os.path.basename(output_folder) + '.ribotyper.long.out')
    if not os.path.exists(ribotyper_long_out):
        cmd = 'ribotyper.pl --skipval -i {cm_library}/modelinfo.txt -f {fasta_input} {output_folder}'.format(
            cm_library=cm_library,
            fasta_input=fasta_input,
            output_folder=output_folder
        )
        print(cmd)
        os.system(cmd)
    f_out = os.path.join(output_folder, 'hits.txt')
    cmd = "cat %s | grep -v '^#' | grep -v MultipleHits | grep PASS | awk -v OFS='\t' '{print $2, $8, $3}' > %s" % (ribotyper_long_out, f_out)
    os.system(cmd)
    return f_out


def symlink_cms(source):
    for cm_file in glob.glob(os.path.join(source, '*.cm')):
        if 'all.cm' not in cm_file:
            target = os.path.join(os.path.abspath(config.CM_LIBRARY), os.path.basename(cm_file))
            if not os.path.exists(target):
                cmd = 'ln -s {} {}'.format(os.path.abspath(cm_file), target)
                os.system(cmd)


@click.group()
def cli():
    pass


@cli.command()
def setup():
    if not os.path.exists(config.CM_LIBRARY):
        os.makedirs(config.CM_LIBRARY)
    crw.setup()
    rfam.setup()
    gtrnadb.setup()
    generate_model_info(cm_library=os.path.join(config.CM_LIBRARY, 'rfam'))


@cli.command()
@click.argument('fasta-input', type=click.Path())
@click.argument('output-folder', type=click.Path())
def draw(fasta_input, output_folder):
    """
    Single entry point for visualising 2D for an RNA sequence.
    Selects a template and runs Traveler using CRW, LSU, or Rfam libraries.
    """
    os.system('mkdir -p %s' % output_folder)

    with open(get_ribotyper_output(fasta_input, output_folder, config.CM_LIBRARY), 'r') as f:
        for line in f.readlines():
            rnacentral_id, model_id, _ = line.split('\t')
            print(line)
            if model_id.count('.') >= 2 or model_id == '5_8S_rRNA':
                crw.visualise_crw(fasta_input, output_folder, rnacentral_id, model_id)
            elif model_id.count('_') == 2:
                ribovision.visualise_lsu(fasta_input, output_folder, rnacentral_id, model_id)
            else:
                rfam.visualise_rfam(fasta_input, output_folder, rnacentral_id, model_id)

    for trna in gtrnadb.classify_trna_sequences(fasta_input, output_folder):
        gtrnadb.generate_2d(trna['domain'], trna['isotype'], trna['id'], trna['start'], trna['end'], fasta_input, output_folder)


@cli.group('gtrnadb')
def gtrnadb_group():
    pass

@gtrnadb_group.command('setup')
def gtrnadb_setup():
    """
    This will copy all the CM files into place so that drawing will not modify
    the data directory.
    """
    gtrnadb.setup()


@gtrnadb_group.command('draw')
@click.option('--test', default=False, is_flag=True, help='Process only the first 10 sequences')
@click.option('--domain', default=False, type=click.STRING, help='Domain (A for Archaea, B for Bacteria, or E for Eukaryotes)')
@click.option('--isotype', default=False, type=click.STRING, help='tRNA isotype, for example Thr')
@click.argument('fasta-input', type=click.Path())
@click.argument('output-folder', type=click.Path())
def gtrnadb_draw(fasta_input, output_folder, domain='', isotype='', test=None):
    """
    Visualise sequences using GtRNAdb templates.
    """
    os.system('mkdir -p %s' % output_folder)

    if domain and isotype:
        gtrnadb.visualise(domain.upper(), isotype.capitalize(), fasta_input, output_folder, test)
    else:
        for trna in gtrnadb.classify_trna_sequences(fasta_input, output_folder):
            gtrnadb.generate_2d(trna['domain'], trna['isotype'], trna['id'], trna['start'], trna['end'], fasta_input, output_folder)


@cli.group('crw')
def crw_group():
    pass


@crw_group.command('draw')
@click.option('--test', default=False, is_flag=True, help='Process only the first 10 sequences')
@click.argument('fasta-input', type=click.Path())
@click.argument('output-folder', type=click.Path())
def rrna_draw(fasta_input, output_folder, test):
    os.system('mkdir -p %s' % output_folder)
    with open(get_ribotyper_output(fasta_input, output_folder, config.CRW_CM_LIBRARY), 'r') as f:
        for line in f.readlines():
            rnacentral_id, model_id, _ = line.split('\t')
            crw.visualise_crw(fasta_input,
                              output_folder,
                              rnacentral_id,
                              model_id)

@cli.group('ribovision')
def ribovision_group():
    """
    Commands dealing with laying out sequences based upon RiboVision models.
    """
    pass


@ribovision_group.command('draw')
@click.argument('fasta-input', type=click.Path())
@click.argument('output-folder', type=click.Path())
def ribovision_draw (fasta_input, output_folder):
    os.system('mkdir -p %s' % output_folder)
    with open(get_ribotyper_output(fasta_input, output_folder, config.RIBOVISION_CM_LIBRARY), 'r') as f:
        for line in f.readlines():
            rnacentral_id, model_id, _ = line.split('\t')
            ribovision.visualise_lsu(fasta_input, output_folder, rnacentral_id, model_id)


@cli.group('rfam')
def rfam_group():
    """
    Commands dealing with laying out sequences based upon Rfam models.
    """
    pass


@rfam_group.command('blacklisted')
def rfam_blacklist():
    """
    Show all blacklisted families. These include rRNA families as well as
    families that do not have any secondary structure.
    """
    for model in sorted(rfam.blacklisted()):
        print(model)


@rfam_group.command('draw')
@click.option('--test', default=False, is_flag=True, help='Process only the first 10 sequences')
@click.argument('rfam_accession', type=click.STRING)
@click.argument('fasta-input', type=click.Path())
@click.argument('output-folder', type=click.Path())
def rfam_draw(rfam_accession, fasta_input, output_folder, test=None):
    """
    Visualise sequences using the Rfam/R-scape consensus structure as template.

    RFAM_ACCESSION - Rfam family to process (RF00001, RF00002 etc)
    """
    print(rfam_accession)
    if rfam_accession == 'all':
        rfam_accs = rfam.get_all_rfam_acc()
    else:
        rfam_accs = [rfam_accession]

    for rfam_acc in rfam_accs:
        if rfam.has_structure(rfam_acc):
            rfam.rscape2traveler(rfam_acc)
            rfam.generate_2d(rfam_acc, output_folder, fasta_input, test)
        else:
            print('{} does not have a conserved secondary structure'.format(rfam_acc))


@rfam_group.command('validate')
@click.argument('rfam_accession', type=click.STRING)
@click.argument('output', type=click.File('w'))
def rfam_validate(rfam_accession, output):
    """
    Check if the given Rfam accession is one that should be drawn. If so it will
    be output to the given file, otherwise it will not.
    """
    if rfam_accession not in rfam.blacklisted():
        output.write(rfam_accession + '\n')


if __name__ == '__main__':
    cli()
