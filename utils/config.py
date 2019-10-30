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
DATA = os.path.join(os.path.dirname(here), 'data')

CRW_CM_LIBRARY = os.path.join(DATA, 'crw-cms')
CRW_PS_LIBRARY = os.path.join(DATA, 'crw-ps')
CRW_BPSEQ_LOCATION = os.path.join(DATA, 'crw-bpseq')
CRW_FASTA_LIBRARY = os.path.join(DATA, 'crw-fasta-no-pseudoknots')

RFAM_DATA = os.path.join(DATA, 'rfam')

RIBOVISION_CM_LIBRARY = os.path.join(DATA, 'ribovision', 'cms')
RIBOVISION_BPSEQ = os.path.join(DATA, 'ribovision', 'bpseq')
RIBOVISION_TRAVELER = os.path.join(DATA, 'ribovision', 'traveler')

CM_LIBRARY = os.path.join(DATA, 'cms')
