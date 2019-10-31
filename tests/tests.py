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

import filecmp
import os
import unittest

from utils import config


EXECUTABLE = os.path.join(config.PROJECT_HOME, 'auto-traveler.py')

# @unittest.skip("")
class TestCovarianceModelDatabase(unittest.TestCase):

    @staticmethod
    def count_lines(filename):
        with open(filename) as f_modelinfo:
            num_lines = sum(1 for line in f_modelinfo)
        return num_lines

    @staticmethod
    def counts_cms(folder):
        return len([name for name in os.listdir(folder)
                    if os.path.isfile(os.path.join(folder,name))
                       and name.endswith('.cm')])

    def verify_cm_database(self, location, count):
        print('Verifying models in {}'.format(location))
        modelinfo = os.path.join(location, 'modelinfo.txt')
        all_cm = os.path.join(location, 'all.cm')
        num_lines = self.count_lines(modelinfo)
        num_cms = self.counts_cms(location)
        self.assertTrue(os.path.exists(modelinfo), 'A required file modelinfo.txt does not exist')
        self.assertTrue(os.path.exists(all_cm), 'A required file all.cm does not exist')
        self.assertEqual(num_lines, num_cms, 'The number of lines in modelinfo.txt file does not match the number of covariance model files')
        self.assertEqual(num_cms, count, 'The number of CMs does not match')

    def test_crw_database(self):
        self.verify_cm_database(config.CRW_CM_LIBRARY, 884)

    def test_ribovision_database(self):
        self.verify_cm_database(config.RIBOVISION_CM_LIBRARY, 18)

    def test_combined_database(self):
        self.verify_cm_database(config.CM_LIBRARY, 3558)

    def test_rfam_database(self):
        pass

# @unittest.skip("")
class TestRibovision(unittest.TestCase):
    fasta_input = os.path.join('examples', 'lsu-small-example.fasta')
    test_results = os.path.join('tests', 'results', 'ribovision')
    precomputed_results = os.path.join('tests', 'examples', 'ribovision')
    cmd = 'python3 {} ribovision draw {} {}'.format(EXECUTABLE, fasta_input, test_results)
    files = [
        'hits.txt',
        'URS000080E357_9606-mHS_LSU_3D.colored.svg',
    ]

    @staticmethod
    def delete_folder(folder):
        os.system('rm -Rf {}'.format(folder))

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(self.cmd)

    def test_examples(self):
        for filename in self.files:
            new_file = os.path.join(self.test_results, filename)
            reference_file = os.path.join(self.precomputed_results, filename)
            self.assertTrue(os.path.exists(new_file))
            self.assertTrue(filecmp.cmp(new_file, reference_file))

    def tearDown(self):
        self.delete_folder(self.test_results)

# @unittest.skip("")
class TestRfam(unittest.TestCase):
    rfam_acc = 'RF00162'
    fasta_input = os.path.join('examples', rfam_acc + '.example.fasta')
    test_results = os.path.join('tests', 'results', 'rfam')
    precomputed_results = os.path.join('tests', 'examples', 'rfam', rfam_acc)
    cmd = 'python3 {} rfam draw {} {} {}'.format(EXECUTABLE, rfam_acc, fasta_input, test_results)
    files = [
        'URS00001D0AD3_224308.colored.svg',
        'URS00002D29F6_224308.colored.svg',
        'URS00002F3927_224308.colored.svg',
        'URS000053CEAC_224308.colored.svg',
        'URS000008638F_224308.colored.svg',
    ]

    @staticmethod
    def delete_folder(folder):
        os.system('rm -Rf {}'.format(folder))

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(self.cmd)

    def test_examples(self):
        for filename in self.files:
            new_file = os.path.join(self.test_results, self.rfam_acc, filename)
            reference_file = os.path.join(self.precomputed_results, filename)
            self.assertTrue(os.path.exists(new_file))
            self.assertTrue(filecmp.cmp(new_file, reference_file))

    def tearDown(self):
        self.delete_folder(self.test_results)

# @unittest.skip("")
class TestCrw(unittest.TestCase):
    label = 'crw'

    fasta_input = os.path.join('examples', label + '-examples.fasta')
    test_results = os.path.join('tests', 'results', label)
    precomputed_results = os.path.join('tests', 'examples', label)
    cmd = 'python3 {} crw draw {} {}'.format(EXECUTABLE, fasta_input, test_results)
    files = [
        'hits.txt',
        'URS00000F9D45_9606-d.5.e.H.sapiens.2.colored.svg',
        'URS000044DFF6_9606-d.16.m.H.sapiens.5.colored.svg',
    ]

    @staticmethod
    def delete_folder(folder):
        os.system('rm -Rf {}'.format(folder))

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(self.cmd)

    def test_examples(self):
        for filename in self.files:
            new_file = os.path.join(self.test_results, filename)
            reference_file = os.path.join(self.precomputed_results, filename)
            self.assertTrue(os.path.exists(new_file))
            self.assertTrue(filecmp.cmp(new_file, reference_file), 'File {} does not match'.format(new_file))

    def tearDown(self):
        self.delete_folder(self.test_results)

# @unittest.skip("")
class TestSingleEntry(unittest.TestCase):
    fasta_input = os.path.join('examples', 'examples.fasta')
    test_results = os.path.join('tests', 'results', 'single-entry')
    precomputed_results = os.path.join('tests', 'examples', 'single-entry')
    cmd = 'python3 {} draw {} {}'.format(EXECUTABLE, fasta_input, test_results)
    files = [
        'hits.txt',
        'URS00000F9D45_9606-d.5.e.H.sapiens.2.colored.svg',
        'URS000044DFF6_9606-d.16.m.H.sapiens.5.colored.svg',
        'URS000053CEAC_224308.colored.svg',
        'URS0000162127_9606.colored.svg',
        'URS000080E357_9606-mHS_LSU_3D.colored.svg',
    ]

    @staticmethod
    def delete_folder(folder):
        os.system('rm -Rf {}'.format(folder))

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(self.cmd)

    def test_examples(self):
        for filename in self.files:
            new_file = os.path.join(self.test_results, filename)
            reference_file = os.path.join(self.precomputed_results, filename)
            self.assertTrue(os.path.exists(new_file))
            self.assertTrue(filecmp.cmp(new_file, reference_file))

    def tearDown(self):
        self.delete_folder(self.test_results)


if __name__ == '__main__':
    unittest.main()
