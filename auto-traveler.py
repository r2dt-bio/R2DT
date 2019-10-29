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

import re
import os

import click

from utils import auto_traveler_rfam as auto_rfam


here = os.path.realpath(os.path.dirname(__file__))
data = os.path.join(here, 'data')

CRW_CM_LIBRARY = os.path.join(data, 'crw-cms')
CRW_PS_LIBRARY = os.path.join(data, 'crw-ps')
CRW_FASTA_LIBRARY = os.path.join(data, 'crw-fasta-no-pseudoknots')
RFAM_DATA = os.path.join(data, 'rfam')

CMS_URL = 'https://www.dropbox.com/s/q5l0s1nj5h4y6e4/cms.tar.gz?dl=0'


def get_ribotyper_output(fasta_input, output_folder, cm_library):
    ribotyper_long_out = os.path.join(output_folder, os.path.basename(output_folder) + '.ribotyper.long.out')
    if not os.path.exists(ribotyper_long_out):
        cmd = 'ribotyper.pl -i {cm_library}/modelinfo.txt -f {fasta_input} {output_folder}'.format(
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


def auto_traveler_rfam(rfam_accession, fasta_input, output_folder, test, rfam_data):
    """
    Visualise sequences using the Rfam/R-scape consensus structure as template.

    RFAM_ACCESSION - Rfam family to process (RF00001, RF00002 etc)
    """
    print(rfam_accession)
    if rfam_accession == 'all':
        rfam_accs = auto_rfam.get_all_rfam_acc()
    else:
        rfam_accs = [rfam_accession]

    for rfam_acc in rfam_accs:
        if auto_rfam.has_structure(rfam_acc):
            auto_rfam.rscape2traveler(rfam_data, rfam_acc)
            auto_rfam.generate_2d(rfam_data, rfam_acc, output_folder, fasta_input, test)
        else:
            print('{} does not have a conserved secondary structure'.format(rfam_acc))


def auto_traveler_lsu(fasta_input, output_folder, test):
    os.system('mkdir -p %s' % output_folder)

    cm_library = os.path.join(data, 'ribovision', 'cms')
    traveler_templates = os.path.join(data, 'ribovision', 'traveler')
    ribovision_bpseq = os.path.join(data, 'ribovision', 'bpseq')

    with open(get_ribotyper_output(fasta_input, output_folder, cm_library), 'r') as f:
        for line in f.readlines():
            rnacentral_id, model_id, _ = line.split('\t')

            cmd = 'esl-sfetch %s %s > temp.fasta' % (fasta_input, rnacentral_id)
            os.system(cmd)

            cmd = "cmalign %s.cm temp.fasta > temp.sto" % os.path.join(cm_library, model_id)
            os.system(cmd)

            cmd = 'esl-alimanip --sindi --outformat pfam temp.sto > temp.stk'
            os.system(cmd)

            cmd = 'ali-pfam-sindi2dot-bracket.pl temp.stk > %s/%s-%s.fasta' % (output_folder, rnacentral_id, model_id)
            os.system(cmd)

            result_base = os.path.join(output_folder, '{rnacentral_id}-{model_id}'.format(
                rnacentral_id=rnacentral_id,
                model_id=model_id,
            ))

            log = result_base + '.log'
            cmd = ('traveler '
                   '--verbose '
                   '--target-structure {result_base}.fasta '
                   '--template-structure --file-format traveler {traveler_templates}/{model_id}.tr {ribovision_bpseq}/{model_id}.fasta '
                   '--all {result_base} > {log}').format(
                       result_base=result_base,
                       model_id=model_id,
                       log=log,
                       traveler_templates=traveler_templates,
                       ribovision_bpseq=ribovision_bpseq
                   )
            print(cmd)
            os.system(cmd)
            os.system('rm temp.fasta temp.sto temp.stk')

            overlaps = 0
            with open(log, 'r') as raw:
                for line in raw:
                    match = re.search(r'Overlaps count: (\d+)', line)
                    if match:
                        if overlaps:
                            print('ERROR: Saw too many overlap counts')
                            break
                        overlaps = int(match.group(1))

            with open(result_base + '.overlaps', 'w') as out:
                out.write(str(overlaps))
                out.write('\n')


def auto_traveler_crw(fasta_input, output_folder, cm_library, ps_library, fasta_library):
    os.system('mkdir -p %s' % output_folder)

    with open(get_ribotyper_output(fasta_input, output_folder, cm_library), 'r') as f:
        for line in f.readlines():
            rnacentral_id, model_id, _ = line.split('\t')
            found = '%s-%s' % (rnacentral_id, model_id)

            cmd = 'esl-sfetch %s %s > temp.fasta' % (fasta_input, rnacentral_id)
            os.system(cmd)

            cmd = "cmalign %s.cm temp.fasta > temp.sto" % os.path.join(cm_library, model_id)
            os.system(cmd)

            cmd = 'esl-alimanip --sindi --outformat pfam temp.sto > temp.stk'
            os.system(cmd)

            cmd = 'ali-pfam-sindi2dot-bracket.pl temp.stk > %s/%s-%s.fasta' % (output_folder, rnacentral_id, model_id)
            os.system(cmd)

            result_base = os.path.join(output_folder, '{rnacentral_id}-{model_id}'.format(
                rnacentral_id=rnacentral_id,
                model_id=model_id,
            ))

            log = result_base + '.log'
            cmd = ('traveler '
                   '--verbose '
                   '--target-structure {result_base}.fasta '
                   '--template-structure {ps_library}/{model_id}.ps {fasta_library}/{model_id}.fasta '
                   '--all {result_base} > {log}').format(
                       result_base=result_base,
                       model_id=model_id,
                       ps_library=ps_library,
                       fasta_library=fasta_library,
                       log=log,
                   )
            print(cmd)
            os.system(cmd)
            os.system('rm temp.fasta temp.sto temp.stk')

            overlaps = 0
            with open(log, 'r') as raw:
                for line in raw:
                    match = re.search(r'Overlaps count: (\d+)', line)
                    if match:
                        if overlaps:
                            print('ERROR: Saw too many overlap counts')
                            break
                        overlaps = int(match.group(1))

            with open(result_base + '.overlaps', 'w') as out:
                out.write(str(overlaps))
                out.write('\n')


@click.group()
def cli():
    pass


@cli.group('rrna')
def rrna_group():
    pass


@rrna_group.command('setup')
@click.option('--cms-url', default=CMS_URL, type=click.STRING)
def rrna_fetch(cms_url=None):
    cms = 'cms.tar.gz'
    cmd = 'wget -O {cms} {url} && tar xf {cms} && rm cms.tar.gz'
    os.system(cmd.format(url=cms_url, cms=cms))

@rrna_group.command('draw')
@click.option('--cm-library', type=click.Path(), default=CRW_CM_LIBRARY)
@click.option('--ps-library', type=click.Path(), default=CRW_PS_LIBRARY)
@click.option('--fasta-library', type=click.Path(), default=CRW_FASTA_LIBRARY)
@click.option('--test', default=False, is_flag=True, help='Process only the first 10 sequences')
@click.argument('fasta-input', type=click.Path())
@click.argument('output-folder', type=click.Path())
def rrna_draw(
    fasta_input,
    output_folder,
    cm_library=None,
    ps_library=None,
    fasta_library=None,
    test=None,
):
    auto_traveler_crw(
        fasta_input,
        output_folder,
        cm_library,
        ps_library,
        fasta_library,
    )

@rrna_group.command('lsu')
@click.argument('fasta-input', type=click.Path())
@click.argument('output-folder', type=click.Path())
def rfam_lsu(fasta_input, output_folder, test=None):
    auto_traveler_lsu(fasta_input, output_folder, test)


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
    for model in sorted(auto_rfam.blacklisted()):
        print(model)


@rfam_group.command('setup')
@click.option('--rfam-data', type=click.Path(), default=RFAM_DATA)
@click.argument('accessions', type=click.STRING, nargs=-1)
def rfam_fetch(accessions, rfam_data=None):
    """
    Fetch data for a given Rfam family. This will be done automatically by the
    pipeline if needed by the drawing step. If given the accession 'all' then
    all Rfam models will be fetched.
    """
    if not accessions:
        accessions = 'all'
    auto_rfam.fetch_data(rfam_data, accessions)


@rfam_group.command('draw')
@click.option('--rfam-data', type=click.Path(), default=RFAM_DATA)
@click.option('--test', default=False, is_flag=True, help='Process only the first 10 sequences')
@click.argument('rfam_accession', type=click.STRING)
@click.argument('fasta-input', type=click.Path())
@click.argument('output-folder', type=click.Path())
def rfam_draw(rfam_accession, fasta_input, output_folder, rfam_data=None, test=None):
    """
    This will draw all sequences in the fasta file using the template from the
    given Rfam family. Files will be produced into the given output folder.
    """
    auto_traveler_rfam(rfam_accession, fasta_input, output_folder, test=test, rfam_data=rfam_data)


@rfam_group.command('validate')
@click.option('--rfam-data', type=click.Path(), default=RFAM_DATA)
@click.argument('rfam_accession', type=click.STRING)
@click.argument('output', type=click.File('w'))
def rfam_validate(rfam_accession, output, rfam_data):
    """
    Check if the given Rfam accession is one that should be drawn. If so it will
    be output to the given file, otherwise it will not.
    """
    if rfam_accession not in auto_rfam.blacklisted():
        output.write(rfam_accession + '\n')


if __name__ == '__main__':
    cli()
