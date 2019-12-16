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
import glob
import re


from . import config
from . import rfam


def generate_model_info(cm_library, rna_type='SSU'):
    print('Processing files in {}'.format(cm_library))

    all_cm = 'all.cm'  # file with all CMs
    all_cm_path = os.path.join(cm_library, all_cm)

    cmd = 'rm -f {all_cm_path} && cat {cm_library}/*.cm > {all_cm_path}'.format(
        cm_library=cm_library,
        all_cm_path=all_cm_path,
    )
    os.system(cmd)

    with open(os.path.join(cm_library, 'modelinfo.txt'), 'w') as f:
        line = '*all*    -    -    %s\n' % all_cm
        f.write(line)
        for cm in glob.glob('%s/*.cm' % cm_library):
            if all_cm in cm:
                continue
            if re.search(r'RF\d{5}', cm):
                with open(cm, 'r') as f_cm:
                    for line in f_cm:
                        if line.startswith('NAME '):
                            model_name = line.strip().split()[-1]
                rfam_acc = os.path.basename(cm).replace('.cm', '')
                if rfam_acc in rfam.BLACKLIST:
                    continue
            else:
                model_name = os.path.basename(cm).replace('.cm', '')
            line = "%s    %s    Bacteria    %s\n" % (model_name, rna_type, os.path.basename(cm))
            f.write(line)
    print('Done')
