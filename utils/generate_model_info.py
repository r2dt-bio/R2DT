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
import shutil
import subprocess as sp

from . import config


def allowed_names(cm_library):
    from . import rfam
    for cm in glob.glob('%s/*.cm' % cm_library):
        model_name = None
        if 'all.cm' in cm:
            continue

        if re.search(r'RF\d{5}', cm):
            rfam_acc = os.path.basename(cm).replace('.cm', '')
            if rfam_acc in rfam.BLACKLIST:
                continue

            with open(cm, 'r') as f_cm:
                for line in f_cm:
                    if line.startswith('NAME '):
                        model_name = line.strip().split()[-1]
        else:
            model_name = os.path.basename(cm).replace('.cm', '')

        if not model_name:
            raise ValueError("Could not find model_name for %s" % cm)
        yield (model_name, cm)


def generate_model_info(cm_library, rna_type='SSU'):
    if not os.path.exists(cm_library):
        raise ValueError("Missing CM directory: " + cm_library)

    print('Processing files in {}'.format(cm_library))

    all_cm_path = os.path.join(cm_library, 'all.cm')
    modelinfo = os.path.join(cm_library, 'modelinfo.txt')

    with open(all_cm_path, 'w') as all_out:
        with open(modelinfo, 'w') as model_out:
            model_out.write('*all*    -    -    all.cm\n')
            for model_name, cm_path in allowed_names(cm_library):
                with open(cm_path, 'r') as cm:
                    shutil.copyfileobj(cm, all_out)
                model_line = "%s    %s    Bacteria    %s\n" % (model_name, rna_type, os.path.basename(cm_path))
                model_out.write(model_line)

    print('Done')
