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


here = os.path.realpath(os.path.dirname(__file__))
PROJECT_HOME = os.path.dirname(here)
DATA = os.path.join(PROJECT_HOME, 'data')

CRW_CM_LIBRARY = os.path.join(DATA, 'cms', 'crw')
CRW_PS_LIBRARY = os.path.join(DATA, 'crw-ps')
CRW_BPSEQ_LOCATION = os.path.join(DATA, 'crw-bpseq')
CRW_FASTA_LIBRARY = os.path.join(DATA, 'crw-fasta-no-pseudoknots')

RFAM_DATA = os.path.join(DATA, 'rfam')

RIBOVISION_LSU_CM_LIBRARY = os.path.join(DATA, 'ribovision-lsu', 'cms')
RIBOVISION_SSU_CM_LIBRARY = os.path.join(DATA, 'ribovision-ssu', 'cms')
RIBOVISION_LSU_BPSEQ = os.path.join(DATA, 'ribovision-lsu', 'bpseq')
RIBOVISION_SSU_BPSEQ = os.path.join(DATA, 'ribovision-ssu', 'bpseq')
RIBOVISION_LSU_TRAVELER = os.path.join(DATA, 'ribovision-lsu', 'traveler')
RIBOVISION_SSU_TRAVELER = os.path.join(DATA, 'ribovision-ssu', 'traveler')

GTRNADB_CM_LIBRARY = os.path.join(DATA, 'cms', 'gtrnadb')
GTRNADB_EUK = os.path.join(DATA, 'gtrnadb', 'eukaryota_isotype_specific')
GTRNADB_BACT = os.path.join(DATA, 'gtrnadb', 'bacteria_isotype_specific')
GTRNADB_ARCH = os.path.join(DATA, 'gtrnadb', 'archaea_isotype_specific')

CM_LIBRARY = os.path.join(DATA, 'cms')
