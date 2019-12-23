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

from . import config
from . import generate_model_info as modelinfo
from . import shared


def setup():
    os.system('rm -Rf {}'.format(config.CRW_CM_LIBRARY))
    cms = os.path.join(config.DATA, 'crw-cms.tar.gz')
    os.system('cd data && tar xf {}'.format(cms))
    modelinfo.generate_model_info(cm_library=config.CRW_CM_LIBRARY)


def visualise_crw(fasta_input, output_folder, rnacentral_id, model_id):

    temp_fasta = tempfile.NamedTemporaryFile()
    temp_sto = tempfile.NamedTemporaryFile()
    temp_stk = tempfile.NamedTemporaryFile()

    cmd = 'esl-sfetch %s %s > %s' % (fasta_input, rnacentral_id, temp_fasta.name)
    os.system(cmd)

    cmd = "cmalign %s.cm %s > %s" % (os.path.join(config.CRW_CM_LIBRARY, model_id), temp_fasta.name, temp_sto.name)
    os.system(cmd)

    cmd = 'esl-alimanip --sindi --outformat pfam {} > {}'.format(temp_sto.name, temp_stk.name)
    os.system(cmd)

    cmd = 'ali-pfam-sindi2dot-bracket.pl %s > %s/%s-%s.fasta' % (temp_stk.name, output_folder, rnacentral_id, model_id)
    os.system(cmd)

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
