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

from . import config


def visualise_lsu(fasta_input, output_folder, rnacentral_id, model_id):
    cmd = 'esl-sfetch %s %s > temp.fasta' % (fasta_input, rnacentral_id)
    os.system(cmd)

    cmd = "cmalign %s.cm temp.fasta > temp.sto" % os.path.join(config.RIBOVISION_CM_LIBRARY, model_id)
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
               traveler_templates=config.RIBOVISION_TRAVELER,
               ribovision_bpseq=config.RIBOVISION_BPSEQ
           )
    print(cmd)
    os.system(cmd)
    os.system('rm -f temp.fasta temp.sto temp.stk')

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
