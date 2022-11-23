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
import subprocess as sp

from . import config
from . import generate_model_info as modelinfo


def setup():
    """Setup CRW CM library."""
    print("Deleting old CRW library")
    os.system(f"rm -Rf {config.CRW_CM_LIBRARY}")
    print("Extracting precomputed CRW archive")
    cmd = ["tar", "xf", "crw-cms.tar.gz"]
    sp.check_output(cmd, cwd=config.DATA)
    cmd = ["mv", "crw-cms", os.path.join(config.CM_LIBRARY, "crw")]
    sp.check_output(cmd, cwd=config.DATA)
    print("Generating CRW modelinfo file")
    modelinfo.generate_model_info(cm_library=config.CRW_CM_LIBRARY)
