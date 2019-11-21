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
import tempfile


# cmfetch /usr/local/lib/tRNAscan-SE/models/TRNAinf-arch-iso arch-Ala


def generate_2d(trnascan_model, fasta_input, output_folder, test):

    destination = '{}/{}'.format(output_folder, trnascan_model)
from . import config


    if not os.path.exists(destination):
        os.makedirs(destination)

    if not os.path.exists(fasta_input + '.ssi'):
        cmd = 'esl-sfetch --index {}'.format(fasta_input)
        os.system(cmd)

    cmd = "grep '>' {} > headers.txt"
    os.system(cmd.format(fasta_input))

    with open('headers.txt', 'r') as f:
        for i, line in enumerate(f):
            if test and i > 10:
                continue
            seq_id = line.split(' ', 1)[0].replace('>', '').strip()
            print(seq_id)
            visualise_trna(fasta_input, destination, seq_id, trnascan_model)
    os.system('rm headers.txt')


def visualise_trna(fasta_input, output_folder, seq_id, trnascan_model):
def get_trnascan_cm(domain, isotype):
    """
    Fetch a domain-specific isotype covariance model as a separate file.
    """
    cm_output = os.path.join(config.GTRNADB_CM_LIBRARY, '{}_{}.cm'.format(domain, isotype))
    if not os.path.exists(cm_output):
        if domain == 'A':
            cm_library = 'TRNAinf-arch-iso'
            cm_name = 'arch-' + isotype
        elif domain == 'B':
            cm_library = 'TRNAinf-bact-iso'
            cm_name = 'bact-' + isotype
        else:
            cm_library = 'TRNAinf-euk-iso'
            cm_name = 'euk-' + isotype
        cmd = 'cmfetch {cm_library} {cm_name} > {cm_output}'.format(
            cm_library=os.path.join('/usr/local/lib/tRNAscan-SE/models', cm_library),
            cm_name=cm_name,
            cm_output=cm_output
        )
        os.system(cmd)
    return cm_output


def get_traveler_template_xml(domain, isotype):
    if domain == 'A':
        return os.path.join(config.GTRNADB_ARCH, 'arch-{}-traveler-template.xml'.format(isotype))
    elif domain == 'B':
        return os.path.join(config.GTRNADB_BACT, 'bact-{}-traveler-template.xml'.format(isotype))
    else:
        return os.path.join(config.GTRNADB_EUK, 'euk-{}-traveler-template.xml'.format(isotype))


def get_traveler_fasta(domain, isotype):
    if domain == 'A':
        return os.path.join(config.GTRNADB_ARCH, 'arch-{}-traveler.fasta'.format(isotype))
    elif domain == 'B':
        return os.path.join(config.GTRNADB_BACT, 'bact-{}-traveler.fasta'.format(isotype))
    else:
        return os.path.join(config.GTRNADB_EUK, 'euk-{}-traveler.fasta'.format(isotype))


    temp_fasta = tempfile.NamedTemporaryFile()
    temp_sto = tempfile.NamedTemporaryFile()
    temp_stk = tempfile.NamedTemporaryFile()

    cmd = 'esl-sfetch %s %s > %s' % (fasta_input, seq_id, temp_fasta.name)
    os.system(cmd)

    cmd = "cmalign {trnascan_cm} {temp_fasta} > {temp_sto}".format(
        trnascan_cm=get_trnascan_cm(domain, isotype),
        temp_fasta=temp_fasta.name,
        temp_sto=temp_sto.name
    )
    os.system(cmd)

    cmd = 'esl-alimanip --sindi --outformat pfam {} > {}'.format(temp_sto.name, temp_stk.name)
    os.system(cmd)

    result_base = os.path.join(output_folder, seq_id.replace('/', '-'))
    input_fasta = os.path.join(output_folder, seq_id + '.fasta')
    cmd = 'ali-pfam-sindi2dot-bracket.pl {} > {}'.format(temp_stk.name, input_fasta)
    os.system(cmd)

    log = result_base + '.log'
    cmd = ('traveler '
           '--verbose '
           '--target-structure {fasta} '
           '--template-structure --file-format traveler {traveler_template_xml} {traveler_fasta} '
           '--all {result_base} '
           '> {log}' ).format(
               fasta=input_fasta,
               result_base=result_base,
               traveler_template_xml=get_traveler_template_xml(domain, isotype),
               traveler_fasta=get_traveler_fasta(domain, isotype),
               log=log
            )
    os.system(cmd)
    print(cmd)

    temp_fasta.close()
    temp_sto.close()
    temp_stk.close()

    cmd = 'rm -f {0}/*.xml {0}/*.ps'.format(output_folder)
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
