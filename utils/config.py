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
import shutil

here = os.path.realpath(os.path.dirname(__file__))
PROJECT_HOME = os.path.dirname(here)
DATA = os.path.join(PROJECT_HOME, "data")


def find_traveler_utils():
    """
    Locate the utils folder of the Traveler installation, which provides
    infernal2mapping.py, enrich_json.py, json2svg.py, and bpseq2fasta.py.
    Candidates are checked in order: the R2DT_TRAVELER_UTILS environment
    variable, the location used in the Docker image, and the installation
    folder of the traveler executable found in PATH. Returns None when no
    candidate contains infernal2mapping.py.
    """
    candidates = []
    if env_path := os.environ.get("R2DT_TRAVELER_UTILS"):
        candidates.append(env_path)
    candidates.append("/rna/traveler/utils")
    if traveler_bin := shutil.which("traveler"):
        traveler_home = os.path.dirname(os.path.dirname(os.path.realpath(traveler_bin)))
        candidates.append(os.path.join(traveler_home, "utils"))
    for candidate in candidates:
        if os.path.isfile(os.path.join(candidate, "infernal2mapping.py")):
            return candidate
    return None


TRAVELER_UTILS = find_traveler_utils()
COLORSCHEME_JSON = os.path.join(here, "colorscheme.json")

CRW_CM_LIBRARY = os.path.join(DATA, "crw")
CRW_PS_LIBRARY = os.path.join(DATA, "crw-ps")
CRW_BPSEQ_LOCATION = os.path.join(DATA, "crw-bpseq")
CRW_FASTA_LIBRARY = os.path.join(DATA, "crw-fasta-no-pseudoknots")

RFAM_DATA = os.path.join(DATA, "rfam")
RFAM_CM_LIBRARY = os.path.join(RFAM_DATA, "cms")

RIBOVISION_LSU = os.path.join(DATA, "ribovision-lsu")
RIBOVISION_SSU = os.path.join(DATA, "ribovision-ssu")
RIBOVISION_LSU_CM_LIBRARY = os.path.join(RIBOVISION_LSU, "cms")
RIBOVISION_SSU_CM_LIBRARY = os.path.join(RIBOVISION_SSU, "cms")
RIBOVISION_LSU_BPSEQ = os.path.join(RIBOVISION_LSU, "bpseq")
RIBOVISION_SSU_BPSEQ = os.path.join(RIBOVISION_SSU, "bpseq")
RIBOVISION_LSU_TRAVELER = os.path.join(RIBOVISION_LSU, "traveler")
RIBOVISION_SSU_TRAVELER = os.path.join(RIBOVISION_SSU, "traveler")

RNASEP = os.path.join(DATA, "rnasep")
RNASEP_CM_LIBRARY = os.path.join(RNASEP, "cms")
RNASEP_BPSEQ = os.path.join(RNASEP, "bpseq")
RNASEP_TRAVELER = os.path.join(RNASEP, "traveler")

TMRNA = os.path.join(DATA, "tmrna")
TMRNA_CM_LIBRARY = os.path.join(TMRNA, "cm")
TMRNA_FASTA_LIBRARY = os.path.join(TMRNA, "fasta")
TMRNA_STO_LIBRARY = os.path.join(TMRNA, "sto")
TMRNA_XML_LIBRARY = os.path.join(TMRNA, "xml")

GTRNADB_CM_LIBRARY = os.path.join(DATA, "gtrnadb", "cms")
GTRNADB_EUK = os.path.join(DATA, "gtrnadb", "eukaryota_isotype_specific")
GTRNADB_BACT = os.path.join(DATA, "gtrnadb", "bacteria_isotype_specific")
GTRNADB_ARCH = os.path.join(DATA, "gtrnadb", "archaea_isotype_specific")
GTRNADB_MITO = os.path.join(DATA, "gtrnadb", "vertebrate_mitochondrial")

LOCAL_DATA = os.path.join(DATA, "local_data")
