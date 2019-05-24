#!/usr/bin/env python

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

from utils.auto_traveler_rfam import get_all_rfam_acc, has_structure, rscape2traveler, generate_2d, echo_blacklist


CM_LIBRARY = '/rna/auto-traveler/data/cms'
CRW_PS_LIBRARY = '/rna/auto-traveler/data/crw-ps'
CRW_FASTA_LIBRARY = '/rna/auto-traveler/data/crw-fasta-no-pseudoknots'


def get_ribotyper_output(fasta_input, output_folder, cm_library):
    ribotyper_long_out = os.path.join(output_folder, os.path.basename(output_folder) + '.ribotyper.long.out')
    if not os.path.exists(ribotyper_long_out):
        cmd = 'ribotyper.pl -i {CM_LIBRARY}/modelinfo.txt -f {fasta_input} {output_folder}'.format(
            CM_LIBRARY=cm_library,
            fasta_input=fasta_input,
            output_folder=output_folder
        )
        print(cmd)
        os.system(cmd)
    f_out = os.path.join(output_folder, 'hits.txt')
    cmd = "cat %s | grep -v '^#' | grep -v MultipleHits | grep PASS | awk -v OFS='\t' '{print $2, $8, $3}' > %s" % (ribotyper_long_out, f_out)
    os.system(cmd)
    return f_out


def auto_traveler_rfam(rfam_accession, fasta_input, output_folder, test):
    """
    Visualise sequences using the Rfam/R-scape consensus structure as template.

    RFAM_ACCESSION - Rfam family to process (RF00001, RF00002 etc)
    """
    print(rfam_accession)
    if rfam_accession == 'all':
        rfam_accs = get_all_rfam_acc()
    else:
        rfam_accs = [rfam_accession]

    for rfam_acc in rfam_accs:
        if has_structure(rfam_acc):
            rscape2traveler(rfam_acc)
            generate_2d(rfam_acc, output_folder, fasta_input, test)
        else:
            print('{} does not have a conserved secondary structure'.format(rfam_acc))


def auto_traveler_ribotyper(fasta_input, output_folder, cm_library, ps_library, fasta_library):
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


@click.command()
@click.option('--fasta-input', type=click.Path(), default='examples/examples.fasta')
@click.option('--output-folder', type=click.Path(), default='temp')
@click.option('--cm-library', type=click.Path(), default=CM_LIBRARY)
@click.option('--ps-library', type=click.Path(), default=CRW_PS_LIBRARY)
@click.option('--fasta-library', type=click.Path(), default=CRW_FASTA_LIBRARY)
@click.option('--rfam-accession', required=False, help='Rfam accession (e.g. RF00162)')
@click.option('--test', default=False, is_flag=True, help='Process only the first 10 sequences')
@click.option('--rfam-blacklist', default=False, is_flag=True, help='Print the list of unsupported Rfam families')
def main(fasta_input, output_folder, test, rfam_blacklist, rfam_accession=None, cm_library=None, ps_library=None, fasta_library=None):
    """
    Generate RNA secondary structure using templates.
    """
    if rfam_blacklist:
        echo_blacklist()
        return

    if rfam_accession:
        auto_traveler_rfam(rfam_accession, fasta_input, output_folder, test)
    else:
        auto_traveler_ribotyper(fasta_input, output_folder, cm_library, ps_library, fasta_library)


if __name__ == '__main__':
    main()
