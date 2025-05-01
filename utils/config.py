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

from pathlib import Path

here = Path(__file__).resolve().parent
PROJECT_HOME = here.parent
DATA = PROJECT_HOME / "data"

CRW_CM_LIBRARY = DATA / "cms" / "crw"
CRW_PS_LIBRARY = DATA / "crw-ps"
CRW_BPSEQ_LOCATION = DATA / "crw-bpseq"
CRW_FASTA_LIBRARY = DATA / "crw-fasta-no-pseudoknots"

RFAM_DATA = DATA / "rfam"

RIBOVISION_LSU = DATA / "ribovision-lsu"
RIBOVISION_SSU = DATA / "ribovision-ssu"
RIBOVISION_LSU_CM_LIBRARY = RIBOVISION_LSU / "cms"
RIBOVISION_SSU_CM_LIBRARY = RIBOVISION_SSU / "cms"
RIBOVISION_LSU_BPSEQ = RIBOVISION_LSU / "bpseq"
RIBOVISION_SSU_BPSEQ = RIBOVISION_SSU / "bpseq"
RIBOVISION_LSU_TRAVELER = RIBOVISION_LSU / "traveler"
RIBOVISION_SSU_TRAVELER = RIBOVISION_SSU / "traveler"

RNASEP = DATA / "rnasep"
RNASEP_CM_LIBRARY = RNASEP / "cms"
RNASEP_BPSEQ = RNASEP / "bpseq"
RNASEP_TRAVELER = RNASEP / "traveler"

TMRNA = DATA / "tmrna"
TMRNA_CM_LIBRARY = TMRNA / "cm"
TMRNA_FASTA_LIBRARY = TMRNA / "fasta"
TMRNA_STO_LIBRARY = TMRNA / "sto"
TMRNA_XML_LIBRARY = TMRNA / "xml"

GTRNADB_CM_LIBRARY = DATA / "gtrnadb" / "cms"
GTRNADB_EUK = DATA / "gtrnadb" / "eukaryota_isotype_specific"
GTRNADB_BACT = DATA / "gtrnadb" / "bacteria_isotype_specific"
GTRNADB_ARCH = DATA / "gtrnadb" / "archaea_isotype_specific"
GTRNADB_MITO = DATA / "gtrnadb" / "vertebrate_mitochondrial"

CM_LIBRARY = DATA / "cms"

LOCAL_DATA = DATA / "local_data"
