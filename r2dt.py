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
import json
import os
import re

import click
from colorhash import ColorHash

from utils import crw, rfam, ribovision, gtrnadb, config, generate_model_info
from utils import generate_model_info as gmi
from utils import list_models as lm


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


def get_seq_ids(input_fasta):
    """
    Get a list of sequence ids from a fasta file.
    """
    seq_ids = set()
    with open(input_fasta, 'r') as f_in:
        for line in f_in:
            if line.startswith('>'):
                match = re.search(r'>(.*?)\s', line)
                if match:
                    seq_ids.add(match.group(1))
    return seq_ids


def get_hits(folder):
    """
    Get a list of sequence ids found in the hits.txt file by ribovore.
    """
    hits = set()
    hits_file = os.path.join(folder, 'hits.txt')
    if not os.path.exists(hits_file):
        return hits
    with open(hits_file, 'r') as f_in:
        for line in f_in:
            hits.add(line.split('\t')[0])
    return hits


def get_subset_fasta(fasta_input, output_filename, seq_ids):
    """
    Extract a fasta file named <output_filename> with sequence ids <seq_ids>
    from <fasta_input>.
    """
    index_filename = output_filename + '.txt'
    with open(index_filename, 'w') as f_out:
        for seq_id in seq_ids:
            f_out.write(seq_id + '\n')
    cmd = 'esl-sfetch -o {} -f {} {}'.format(output_filename, fasta_input, index_filename)
    os.system(cmd)
    os.system('esl-sfetch --index ' + output_filename)


@cli.command()
@click.argument('fasta-input', type=click.Path())
@click.argument('output-folder', type=click.Path())
@click.option('--force_template', type=click.STRING, default=None, help='Force sequences into a specific template')
@click.pass_context
def draw(ctx, fasta_input, output_folder, force_template):
    """
    Single entry point for visualising 2D for an RNA sequence.
    Selects a template and runs Traveler using CRW, LSU, or Rfam libraries.
    """
    all_seq_ids = get_seq_ids(fasta_input)

    if force_template:
        for seq_id in all_seq_ids:
            force_draw(force_template, fasta_input, output_folder, seq_id)
        return

    os.system('mkdir -p %s' % output_folder)
    crw_output = os.path.join(output_folder, 'crw')
    ribovision_ssu_output = os.path.join(output_folder, 'ribovision-ssu')
    ribovision_lsu_output = os.path.join(output_folder, 'ribovision-lsu')
    rfam_output = os.path.join(output_folder, 'rfam')
    gtrnadb_output = os.path.join(output_folder, 'gtrnadb')
    rfam_trna_output = os.path.join(output_folder, 'RF00005')
    rnasep_output = os.path.join(output_folder, 'rnasep')

    hits = set()
    subset_fasta = os.path.join(output_folder, 'subset.fasta')
    os.system('cp {} {}'.format(fasta_input, subset_fasta))
    os.system('esl-sfetch --index ' + subset_fasta)

    # Rfam
    print('Analysing {} sequences with Rfam'.format(len(all_seq_ids)))
    with open(get_ribotyper_output(fasta_input, rfam_output, os.path.join(config.CM_LIBRARY, 'rfam')), 'r') as f:
        for line in f.readlines():
            rnacentral_id, model_id, _ = line.split('\t')
            rfam.visualise_rfam(fasta_input, rfam_output, rnacentral_id, model_id)

    # RiboVision SSU
    hits = hits.union(get_hits(rfam_output))
    subset = all_seq_ids.difference(hits)
    if subset:
        get_subset_fasta(fasta_input, subset_fasta, subset)
        print('Analysing {} sequences with RiboVision SSU'.format(len(subset)))
        ctx.invoke(ribovision_draw_ssu, fasta_input=subset_fasta, output_folder=ribovision_ssu_output)

    # CRW
    hits = hits.union(get_hits(ribovision_ssu_output))
    subset = all_seq_ids.difference(hits)
    if subset:
        get_subset_fasta(fasta_input, subset_fasta, subset)
        print('Analysing {} sequences with CRW'.format(len(subset)))
        ctx.invoke(rrna_draw, fasta_input=subset_fasta, output_folder=crw_output, test=False)

    # RiboVision LSU
    hits = hits.union(get_hits(crw_output))
    subset = all_seq_ids.difference(hits)
    if subset:
        get_subset_fasta(fasta_input, subset_fasta, subset)
        print('Analysing {} sequences with RiboVision LSU'.format(len(subset)))
        ctx.invoke(ribovision_draw_lsu, fasta_input=subset_fasta, output_folder=ribovision_lsu_output)

    # RNAse P
    hits = hits.union(get_hits(ribovision_lsu_output))
    subset = all_seq_ids.difference(hits)
    if subset:
        get_subset_fasta(fasta_input, subset_fasta, subset)
        print('Analysing {} sequences with RNAse P models'.format(len(subset)))
        ctx.invoke(rnasep_draw, fasta_input=subset_fasta, output_folder=rnasep_output)

    # GtRNAdb
    hits = hits.union(get_hits(rnasep_output))
    subset = all_seq_ids.difference(hits)
    if subset:
        get_subset_fasta(fasta_input, subset_fasta, subset)
        print('Analysing {} sequences with GtRNAdb'.format(len(subset)))
        for trna in gtrnadb.classify_trna_sequences(subset_fasta, gtrnadb_output):
            gtrnadb.generate_2d(trna['domain'], trna['isotype'], trna['id'], trna['start'], trna['end'], fasta_input, output_folder + '/gtrnadb')

    # Rfam tRNA
    hits = hits.union(get_hits(gtrnadb_output))
    subset = all_seq_ids.difference(hits)
    if subset:
        get_subset_fasta(fasta_input, subset_fasta, subset)
        print('Analysing {} sequences with Rfam tRNA'.format(len(subset)))
        trna_ids = rfam.cmsearch_nohmm_mode(subset_fasta, output_folder, 'RF00005')
        if trna_ids:
            get_subset_fasta(fasta_input, subset_fasta, trna_ids)
            rfam.generate_2d('RF00005', output_folder, subset_fasta, False)

    # move svg files to the final location
    result_folders = [crw_output, ribovision_ssu_output, ribovision_lsu_output, rfam_output, gtrnadb_output, rfam_trna_output, rnasep_output]
    for folder in result_folders:
        organise_results(folder, output_folder)
    organise_metadata(output_folder, result_folders)

    # clean up
    os.system('rm {}/subset*'.format(output_folder))


def organise_results(results_folder, output_folder):
    destination = os.path.join(output_folder, 'results')
    svg_folder = os.path.join(destination, 'svg')
    thumbnail_folder = os.path.join(destination, 'thumbnail')
    fasta_folder = os.path.join(destination, 'fasta')
    for folder in [destination, svg_folder, thumbnail_folder, fasta_folder]:
        os.system('mkdir -p {}'.format(folder))

    svgs = glob.glob(os.path.join(results_folder, '*.colored.svg'))
    if len(svgs):
        for svg in svgs:
            with open(svg, 'r') as f_svg:
                thumbnail = generate_thumbnail(f_svg.read(), svg)
                with open(svg.replace('.colored.', '.thumbnail.'), 'w') as f_thumbnail:
                    f_thumbnail.write(thumbnail)
        os.system('mv {0}/*.colored.svg {1}'.format(results_folder, svg_folder))
        os.system('mv {0}/*.thumbnail.svg {1}'.format(results_folder, thumbnail_folder))
        os.system('mv {0}/*.fasta {1}'.format(results_folder, fasta_folder))


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


@cli.group('rnasep')
def rnasep_group():
    pass

@rnasep_group.command('draw')
@click.argument('fasta-input', type=click.Path())
@click.argument('output-folder', type=click.Path())
def rnasep_draw(fasta_input, output_folder):
    os.system('mkdir -p %s' % output_folder)
    with open(get_ribotyper_output(fasta_input, output_folder, config.RNASEP_CM_LIBRARY), 'r') as f:
        for line in f.readlines():
            rnacentral_id, model_id, _ = line.split('\t')
            ribovision.visualise('rnasep', fasta_input, output_folder, rnacentral_id, model_id)


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


@ribovision_group.command('draw_lsu')
@click.argument('fasta-input', type=click.Path())
@click.argument('output-folder', type=click.Path())
def ribovision_draw_lsu(fasta_input, output_folder):
    os.system('mkdir -p %s' % output_folder)
    with open(get_ribotyper_output(fasta_input, output_folder, config.RIBOVISION_LSU_CM_LIBRARY), 'r') as f:
        for line in f.readlines():
            rnacentral_id, model_id, _ = line.split('\t')
            ribovision.visualise('lsu', fasta_input, output_folder, rnacentral_id, model_id)


@ribovision_group.command('draw_ssu')
@click.argument('fasta-input', type=click.Path())
@click.argument('output-folder', type=click.Path())
def ribovision_draw_ssu(fasta_input, output_folder):
    os.system('mkdir -p %s' % output_folder)
    with open(get_ribotyper_output(fasta_input, output_folder, config.RIBOVISION_SSU_CM_LIBRARY), 'r') as f:
        for line in f.readlines():
            rnacentral_id, model_id, _ = line.split('\t')
            ribovision.visualise('ssu', fasta_input, output_folder, rnacentral_id, model_id)


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


def generate_thumbnail(image, description):
    move_to_start_position = None
    color = ColorHash(description).hex
    points = []
    for i, line in enumerate(image.split('\n')):
        if 'width' in line and not 'stroke-width' in line:
            width = re.findall(r'width="(\d+(\.\d+)?)"', line)
        if 'height' in line:
            height = re.findall(r'height="(\d+(\.\d+)?)"', line)
        for nt in re.finditer('<text x="(\d+)(\.\d+)?" y="(\d+)(\.\d+)?".*?</text>', line):
            if 'numbering-label' in nt.group(0):
                continue
            if not move_to_start_position:
                move_to_start_position = 'M{} {} '.format(nt.group(1), nt.group(3))
            points.append('L{} {}'.format(nt.group(1), nt.group(3)))
    if len(points) < 200:
        stroke_width = '3'
    elif len(points) < 500:
        stroke_width = '4'
    elif len(points) < 3000:
        stroke_width = '4'
    else:
        stroke_width = '2'
    thumbnail = '<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}"><path style="stroke:{};stroke-width:{}px;fill:none;" d="'.format(width[0][0], height[0][0], color, stroke_width)
    thumbnail += move_to_start_position
    thumbnail += ' '.join(points)
    thumbnail += '"/></svg>'
    return thumbnail


def organise_metadata(output_folder, result_folders):
    """
    Aggregate hits.txt files from all subfolders.
    """
    tsv_folder = os.path.join(output_folder, 'results', 'tsv')
    os.system('mkdir -p {}'.format(tsv_folder))
    with open(os.path.join(tsv_folder, 'metadata.tsv'), 'w') as f_out:
        for folder in result_folders:
            hits = os.path.join(folder, 'hits.txt')
            if not os.path.exists(hits):
                continue
            with open(hits, 'r') as f_hits:
                for line in f_hits.readlines():
                    if 'gtrnadb' in folder:
                        line = line.replace('PASS', 'GtRNAdb')
                    elif 'crw' in folder:
                        line = line.replace('PASS', 'CRW')
                    elif 'rfam' in folder or 'RF00005' in folder:
                        line = line.replace('PASS', 'Rfam')
                    elif 'ribovision-lsu' in folder or 'ribovision-ssu' in folder:
                        line = line.replace('PASS', 'RiboVision')
                    elif 'rnasep' in folder:
                        line = line.replace('PASS', 'RNAse P Database')
                    f_out.write(line)


@cli.command()
@click.argument('cm_library', type=click.Path())
def generatemodelinfo(cm_library):
    gmi.generate_model_info(cm_library)


def force_draw(model_id, fasta_input, output_folder, seq_id):
    model_type = lm.get_model_type(model_id)
    if not model_type:
        print('Error: Model not found. Please check model_id')
        return
    print('Visualising sequence {} using the {} model from {}'.format(seq_id, model_id, model_type))
    os.system('esl-sfetch --index %s' % fasta_input)

    output = os.path.join(output_folder, model_type.replace('_', '-'))

    if model_type == 'rfam':
        rfam.visualise_rfam(fasta_input, output, seq_id, model_id)
    elif model_type == 'ribovision_ssu':
        ribovision.visualise('ssu', fasta_input, output, seq_id, model_id)
    elif model_type == 'ribovision_lsu':
        ribovision.visualise('lsu', fasta_input, output, seq_id, model_id)
    elif model_type == 'rnasep':
        ribovision.visualise('rnasep', fasta_input, output, seq_id, model_id)
    elif model_type == 'crw':
        crw.visualise_crw(fasta_input, output, seq_id, model_id)
    elif model_type == 'gtrnadb':
        domain, isotype = model_id.split('_')
        test = False
        gtrnadb.visualise(domain, isotype, fasta_input, output, test)
    # organise results into folders
    organise_results(output, output_folder)
    metadata_folder = os.path.join(output_folder, 'results', 'tsv')
    if not os.path.exists(metadata_folder):
        os.makedirs(metadata_folder)
    label_mapping = {
        'crw': 'CRW',
        'gtrnadb': 'GtRNAdb',
        'rfam': 'Rfam',
        'ribovision_ssu': 'RiboVision',
        'ribovision_lsu': 'RiboVision',
        'rnasep': 'RNAse P database',
    }
    with open(os.path.join(metadata_folder, 'metadata.tsv'), 'a') as f_out:
        line = '{}\t{}\t{}\n'.format(seq_id, model_id, label_mapping[model_type])
        f_out.write(line)


@cli.command()
def list_models():
    data = lm.list_models()
    for item in data:
        print(item['description'])
    with open(os.path.join(config.DATA, 'models.json'), 'w') as models_file:
        json.dump(data, models_file)

if __name__ == '__main__':
    cli()
