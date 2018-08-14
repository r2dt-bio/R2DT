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


CM_LIBRARY = '/rna/auto-traveler/data/crw-cm'


all_cm = 'all.cm'  # file with all CMs
all_cm_path = os.path.join(CM_LIBRARY, all_cm)

cmd = 'rm {all_cm_path} && cat {CM_LIBRARY}/*.cm > {all_cm_path}'.format(CM_LIBRARY=CM_LIBRARY, all_cm_path=all_cm_path)
os.system(cmd)

with open(os.path.join(CM_LIBRARY, 'modelinfo.txt'), 'w') as f:
    line = '*all*    -    -    %s\n' % all_cm
    f.write(line)
    for cm in glob.glob('%s/*.cm' % CM_LIBRARY):
        if all_cm in cm:
            continue
        model_name = os.path.basename(cm).replace('.cm', '')
        line = "%s    SSU    Bacteria    %s\n" % (model_name, os.path.basename(cm))
        f.write(line)
print 'Done'
