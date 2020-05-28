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
import subprocess as sp

from . import config
from . import generate_model_info as modelinfo
from . import shared


def setup():
    print('Deleting old CRW library')
    os.system('rm -Rf {}'.format(config.CRW_CM_LIBRARY))
    print('Extracting precomputed CRW archive')
    cmd = ['tar', 'xf', 'crw-cms.tar.gz']
    sp.check_output(cmd, cwd=config.DATA)
    cmd = ['mv', 'crw-cms', os.path.join(config.CM_LIBRARY, 'crw')]
    sp.check_output(cmd, cwd=config.DATA)
    print('Generating CRW modelinfo file')
    modelinfo.generate_model_info(cm_library=config.CRW_CM_LIBRARY)


def visualise_crw(fasta_input, output_folder, rnacentral_id, model_id):

    temp_fasta = tempfile.NamedTemporaryFile()
    temp_sto = tempfile.NamedTemporaryFile()
    temp_stk = tempfile.NamedTemporaryFile()

    cmd = 'esl-sfetch %s %s > %s' % (fasta_input, rnacentral_id, temp_fasta.name)
    result = os.system(cmd)
    if result:
        raise ValueError("Failed esl-sfetch for: %s" % rnacentral_id)

    model_path = os.path.join(config.CRW_CM_LIBRARY, model_id)
    if not os.path.exists(model_path + '.cm'):
        print('Model %s does not exist' % model_path)
        return
    cm_options = ['', '--cyk --notrunc --noprob --nonbanded --small']
    for options in cm_options:
        cmd = "cmalign %s %s.cm %s > %s" % (options, model_path, temp_fasta.name, temp_sto.name)
        result = os.system(cmd)
        if not result:
            break
    else:
        print("Failed cmalign of %s to %s" % (rnacentral_id, model_id))
        return

    cmd = 'esl-alimanip --sindi --outformat pfam {} > {}'.format(temp_sto.name, temp_stk.name)
    result = os.system(cmd)
    if result:
        raise ValueError("Failed esl-alimanip for %s %s" % (rnacentral_id, model_id))

    cmd = 'ali-pfam-sindi2dot-bracket.pl %s > %s/%s-%s.fasta' % (temp_stk.name, output_folder, rnacentral_id, model_id)
    result = os.system(cmd)
    if result:
        raise ValueError("Failed esl-pfam-sindi2dot-bracket for %s %s" % (rnacentral_id, model_id))

    result_base = os.path.join(output_folder, '{rnacentral_id}-{model_id}'.format(
        rnacentral_id=rnacentral_id,
        model_id=model_id,
    ))

    shared.remove_large_insertions(result_base + '.fasta')

    log = result_base + '.log'
    cmd = ('traveler '
           '--verbose '
           '--target-structure {result_base}.fasta '
           '--template-structure {ps_library}/{model_id}.ps {fasta_library}/{model_id}.fasta '
           '--all {result_base} > {log}').format(
               result_base=result_base,
               model_id=model_id,
               ps_library=config.CRW_PS_LIBRARY,
               fasta_library=config.CRW_FASTA_LIBRARY,
               log=log,
           )
    print(cmd)
    os.system(cmd)

    temp_fasta.close()
    temp_sto.close()
    temp_stk.close()

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
