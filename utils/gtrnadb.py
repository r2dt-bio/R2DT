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
    temp_fasta = tempfile.NamedTemporaryFile()
    temp_sto = tempfile.NamedTemporaryFile()
    temp_stk = tempfile.NamedTemporaryFile()

    cmd = 'esl-sfetch %s %s > %s' % (fasta_input, seq_id, temp_fasta.name)
    os.system(cmd)

    # /rna/tRNAscan-SE-2.0/lib/models/TRNAinf-euk.cm
    trnascan_cm = '/rna/tRNAscan-SE-2.0/lib/models/{}.cm'.format(trnascan_model)
    cmd = "cmalign {trnascan_cm} {temp_fasta} > {temp_sto}".format(
        trnascan_cm=trnascan_cm,
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
               traveler_template_xml='/rna/auto-traveler/data/gtrnadb/eukaryota_isotype_specific/euk-Thr-traveler-template.xml',
               traveler_fasta='/rna/auto-traveler/data/gtrnadb/eukaryota_isotype_specific/euk-Thr-traveler.fasta',
               log=log
            )
    os.system(cmd)

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
