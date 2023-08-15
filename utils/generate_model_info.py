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


import glob
import os
import re
import shutil

from . import rfam


def allowed_names(cm_library):
    """Return all model names and CM paths
    for a given template library."""
    for cm_file in glob.glob(os.path.join(cm_library, "*.cm")):
        model_name = None
        if "all.cm" in cm_file:
            continue
        if re.search(r"RF\d{5}", cm_file):
            rfam_acc = os.path.basename(cm_file).replace(".cm", "")
            if rfam_acc in rfam.BLACKLIST:
                continue
            with open(cm_file, "r", encoding="utf-8") as f_cm:
                for line in f_cm:
                    if line.startswith("NAME "):
                        model_name = line.strip().split()[-1]
        else:
            model_name = os.path.basename(cm_file).replace(".cm", "")
        if not model_name:
            raise ValueError(f"Could not find model_name for {cm_file}")
        yield (model_name, cm_file)


def generate_model_info(cm_library, rna_type="SSU"):
    """Generate a model info file listing all covariance models
    for a given template library."""
    if not os.path.exists(cm_library):
        raise ValueError("Missing CM directory: " + cm_library)
    print(f"Processing files in {cm_library}")

    all_cm_path = os.path.join(cm_library, "all.cm")
    modelinfo = os.path.join(cm_library, "modelinfo.txt")

    with open(all_cm_path, "w", encoding="utf-8") as all_out:
        with open(modelinfo, "w", encoding="utf-8") as model_out:
            model_out.write("*all*    -    -    all.cm\n")
            for model_name, cm_path in allowed_names(cm_library):
                with open(cm_path, "r", encoding="utf-8") as cm_file:
                    shutil.copyfileobj(cm_file, all_out)
                model_line = "    ".join(
                    [model_name, rna_type, "Bacteria", os.path.basename(cm_path)]
                )
                model_out.write(f"{model_line}\n")
    print("Done")
