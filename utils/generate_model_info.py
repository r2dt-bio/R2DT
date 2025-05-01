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


import re
import shutil
from pathlib import Path

from . import rfam


def allowed_names(cm_library):
    """Return all model names and CM paths
    for a given template library."""
    cm_library = Path(cm_library)
    for cm_file in cm_library.glob("*.cm"):
        model_name = None
        if cm_file.name == "all.cm":
            continue
        if re.search(r"RF\d{5}", str(cm_file)):
            rfam_acc = cm_file.stem  # removes .cm
            if rfam_acc in rfam.BLACKLIST:
                continue
            with cm_file.open("r", encoding="utf-8") as f_cm:
                for line in f_cm:
                    if line.startswith("NAME "):
                        model_name = line.strip().split()[-1]
        else:
            model_name = cm_file.stem
        if not model_name:
            raise ValueError(f"Could not find model_name for {cm_file}")
        yield (model_name, cm_file)


def generate_model_info(cm_library, rna_type="SSU"):
    """Generate a model info file listing all covariance models
    for a given template library."""
    cm_library = Path(cm_library)
    if not cm_library.exists():
        raise ValueError("Missing CM directory: " + str(cm_library))
    print(f"Processing files in {cm_library}")

    all_cm_path = cm_library / "all.cm"
    modelinfo = cm_library / "modelinfo.txt"

    with all_cm_path.open("w", encoding="utf-8") as all_out:
        with modelinfo.open("w", encoding="utf-8") as model_out:
            model_out.write("*all*    -    -    all.cm\n")
            for model_name, cm_path in allowed_names(cm_library):
                with cm_path.open("r", encoding="utf-8") as cm_file:
                    shutil.copyfileobj(cm_file, all_out)
                model_line = "    ".join(
                    [model_name, rna_type, "Bacteria", cm_path.name]
                )
                model_out.write(f"{model_line}\n")
    print("Done")
