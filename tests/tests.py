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

from utils import config, rfam


class R2dtTestCase(unittest.TestCase):
    """Base class for R2DT tests."""

    test_results = "test_folder"
    test_results_subfolder = ""
    files = []
    precomputed_results = ""

    @staticmethod
    def delete_folder(folder):
        """Delete test folder"""
        os.system(f"rm -Rf {folder}")

    def tearDown(self):
        if os.environ.get("R2DT_KEEP_TEST_RESULTS", False) == "1":
            print(f"Test results can be found in {self.test_results}")
        else:
            self.delete_folder(self.test_results)

    @staticmethod
    def create_webpage(filename, before, after):
        """Create an HTML file comparing the reference SVG with a new one."""
        with open(filename, "w", encoding="utf-8") as f_html:
            f_html.write("<html><body>")
            with open(before, "r", encoding="utf-8") as f_before:
                svg = f_before.read()
            f_html.write(svg)
            with open(after, "r", encoding="utf-8") as f_after:
                svg = f_after.read()
            f_html.write(svg)
            f_html.write("</body></html>")

    def check_examples(self):
        """Check that files exist and are identical to examples."""
        count = 0
        html_files = []
        for loop_id, filename in enumerate(self.files):
            new_file = os.path.join(
                self.test_results, self.test_results_subfolder, filename
            )
            reference_file = os.path.join(self.precomputed_results, filename)
            self.assertTrue(os.path.exists(new_file))
            is_identical = filecmp.cmp(new_file, reference_file)
            if not is_identical:
                filename = os.path.join(
                    "tests", f"{self.__class__.__name__}_{loop_id}.html"
                )
                self.create_webpage(filename, reference_file, new_file)
                count += 1
                html_files.append(filename)
        self.assertEqual(
            count, 0, f"Found {count} out of {len(self.files)} non-identical SVG files"
        )
        if html_files:
            print(f"Please inspect {', '.join(html_files)}")


class TestCovarianceModelDatabase(unittest.TestCase):
    """Check that the covariance model and template files exist."""

    @staticmethod
    def count_lines(filename):
        """Return the number of lines in a modelinfo file."""
        with open(filename, encoding="utf-8") as f_modelinfo:
            num_lines = sum(1 for line in f_modelinfo)
        return num_lines

    @staticmethod
    def counts_cms(folder):
        """Return the number of .cm files in a folder."""
        return len(
            [
                name
                for name in os.listdir(folder)
                if os.path.isfile(os.path.join(folder, name)) and name.endswith(".cm")
            ]
        )

    def verify_cm_database(self, location, count):
        """Check that the required files exist."""
        print(f"Verifying models in {location}")
        modelinfo = os.path.join(location, "modelinfo.txt")
        all_cm = os.path.join(location, "all.cm")
        num_lines = self.count_lines(modelinfo)
        num_cms = self.counts_cms(location)
        self.assertTrue(
            os.path.exists(modelinfo), "A required file modelinfo.txt does not exist"
        )
        self.assertTrue(os.path.exists(all_cm), "A required file all.cm does not exist")
        self.assertEqual(
            num_lines,
            num_cms,
            "The number of lines in modelinfo.txt does not match the number of covariance models",
        )
        self.assertEqual(num_cms, count, "The number of CMs does not match")

    def test_crw_database(self):
        """Check CRW covariance models."""
        self.verify_cm_database(config.CRW_CM_LIBRARY, 884)

    def test_ribovision_lsu_database(self):
        """Check RiboVision LSU covariance models."""
        self.verify_cm_database(config.RIBOVISION_LSU_CM_LIBRARY, 22)

    def test_ribovision_ssu_database(self):
        """Check RiboVision SSU covariance models."""
        self.verify_cm_database(config.RIBOVISION_SSU_CM_LIBRARY, 11)

    def test_rnasep_cm_database(self):
        """Check RNAse P covariance models."""
        self.verify_cm_database(config.RNASEP_CM_LIBRARY, 24)

    def test_rfam_database(self):
        """
        Check Rfam covariance models and templates.
        As Rfam templates are generated automatically with each
        Rfam release, it is necessary to perform more in depth
        checks than for other template sources that do not change.
        """
        for rfam_acc in rfam.get_all_rfam_acc():
            if rfam_acc in rfam.blacklisted():
                continue
            template = rfam.get_traveler_template_xml(rfam_acc)
            self.assertTrue(os.path.exists(template), f"{template} not found")
            fasta = rfam.get_traveler_fasta(rfam_acc)
            self.assertTrue(os.path.exists(fasta), f"{fasta} not found")
            cm_file = rfam.get_rfam_cm(rfam_acc)
            self.assertTrue(os.path.exists(cm_file), f"{cm_file} not found")


class TestRibovisionLSU(R2dtTestCase):
    """Check that the RiboVision LSU templates work."""

    fasta_input = os.path.join("examples", "lsu-small-example.fasta")
    test_results = os.path.join("tests", "results", "ribovision")
    precomputed_results = os.path.join("tests", "examples", "ribovision")
    cmd = f"r2dt.py ribovision draw_lsu {fasta_input} {test_results}"
    files = [
        "hits.txt",
        "URS000080E357_9606-mHS_LSU_3D.colored.svg",
    ]

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(self.cmd)

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestRibovisionSSU(R2dtTestCase):
    """Check that the RiboVision SSU templates work."""

    fasta_input = os.path.join("examples", "ribovision-ssu-examples.fasta")
    test_results = os.path.join("tests", "results", "ribovision-ssu")
    precomputed_results = os.path.join("tests", "examples", "ribovision-ssu")
    cmd = f"r2dt.py ribovision draw_ssu {fasta_input} {test_results}"
    files = [
        "hits.txt",
        "URS00002A2E83_10090-HS_SSU_3D.colored.svg",
    ]

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(self.cmd)

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestRfamAccession(R2dtTestCase):
    """Test Rfam visualisation when specifying an Rfam accession."""

    rfam_acc = "RF00162"
    fasta_input = os.path.join("examples", rfam_acc + ".example.fasta")
    test_results = os.path.join("tests", "results", "rfam")
    test_results_subfolder = rfam_acc
    precomputed_results = os.path.join("tests", "examples", "rfam", rfam_acc)
    cmd = f"r2dt.py rfam draw {rfam_acc} {fasta_input} {test_results}"
    files = [
        "URS00001D0AD3_224308-RF00162.colored.svg",
        "URS00002D29F6_224308-RF00162.colored.svg",
        "URS00002F3927_224308-RF00162.colored.svg",
        "URS000053CEAC_224308-RF00162.colored.svg",
        "URS000008638F_224308-RF00162.colored.svg",
    ]

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(self.cmd)

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestRfam(R2dtTestCase):
    """
    Sequences from different Rfam families with pseudoknots.
    """

    fasta_input = os.path.join("examples", "rfam.fasta")
    test_results = os.path.join("tests", "results", "rfam", "combined")
    test_results_subfolder = os.path.join("results", "svg")
    precomputed_results = os.path.join("tests", "examples", "rfam", "combined")
    cmd = f"r2dt.py draw {fasta_input} {test_results}"
    files = [
        "URS00021ED9B3_2697049-RF00507.colored.svg",
        "URS000080E2F0_93929-RF01734.colored.svg",
        "URS0000868535_32630-RF01750.colored.svg",
    ]

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(self.cmd)

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestCrw(R2dtTestCase):
    """Check that the CRW templates work."""

    label = "crw"
    fasta_input = os.path.join("examples", label + "-examples.fasta")
    test_results = os.path.join("tests", "results", label)
    precomputed_results = os.path.join("tests", "examples", label)
    cmd = f"r2dt.py crw draw {fasta_input} {test_results}"
    files = [
        "hits.txt",
        "URS00000F9D45_9606-d.5.e.H.sapiens.2.colored.svg",
        "URS000044DFF6_9606-d.16.m.H.sapiens.5.colored.svg",
        "URS000001AE2D_4932-d.16.e.S.cerevisiae.colored.svg",
    ]

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(self.cmd)

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestSingleEntry(R2dtTestCase):
    """Check a broad range of RNA types at once."""

    fasta_input = os.path.join("examples", "examples.fasta")
    test_results = os.path.join("tests", "results", "single-entry")
    test_results_subfolder = os.path.join("results", "svg")
    precomputed_results = os.path.join("tests", "examples", "single-entry")
    cmd = f"r2dt.py draw {fasta_input} {test_results}"
    files = [
        "URS00000F9D45_9606-d.5.e.H.sapiens.2.colored.svg",
        "URS000044DFF6_9606-d.16.m.H.sapiens.5.colored.svg",
        "URS000053CEAC_224308-RF00162.colored.svg",
        "URS0000162127_9606-RF00003.colored.svg",
        "URS000080E357_9606-mHS_LSU_3D.colored.svg",
        "URS0000023412_9606-E_Thr.colored.svg",
        "URS00000012EC-M_Ile.colored.svg",
        "URS0000664B0C_4896-RNAseP_e_S_pombe_JB.colored.svg",
    ]

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(self.cmd)

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()

    def test_json_files(self):
        """Check that files exist and are identical to examples."""
        for filename in self.files:
            filename = filename.replace("svg", "json")
            json_file = os.path.join(self.test_results, "results", "json", filename)
            self.assertTrue(
                os.path.exists(json_file),
                f"Json file {json_file} does not exist",
            )


class TestGtrnadbDomainIsotype(R2dtTestCase):
    """Check tRNA visualisation when specifying both domain and isotype."""

    trnascan_model = "E_Thr"
    fasta_input = os.path.join("examples", f"gtrnadb.{trnascan_model}.fasta")
    test_results = os.path.join("tests", "results", "gtrnadb")
    test_results_subfolder = trnascan_model
    precomputed_results = os.path.join("tests", "examples", "gtrnadb", trnascan_model)
    cmd = f"r2dt.py gtrnadb draw {fasta_input} {test_results} --domain E --isotype Thr"
    files = [
        "URS0000023412_9606-E_Thr.colored.svg",
        "URS000021550A_9606-E_Thr.colored.svg",
        "URS00000A1A88_9606-E_Thr.colored.svg",
        "URS00000F30A4_9606-E_Thr.colored.svg",
        "URS00001D9AFB_9606-E_Thr.colored.svg",
    ]

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(self.cmd)

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestGtrnadbMitoVert(R2dtTestCase):
    """Test mitochondrial vertebrate tRNA templates."""

    fasta_input = os.path.join("examples", "gtrnadb-mito-vert.fasta")
    test_results = os.path.join("tests", "results", "gtrnadb", "mito-vert")
    precomputed_results = os.path.join("tests", "examples", "gtrnadb", "mito-vert")
    cmd = f"r2dt.py gtrnadb draw {fasta_input} {test_results}"
    files = [
        "URS000061A10B_9606-M_LeuTAA.colored.svg",
        "URS000054F2AC_109923-M_LeuTAG.colored.svg",
        "URS0000333A94_392897-M_SerTGA.colored.svg",
        "URS0000043FFB_392897-M_SerGCT.colored.svg",
        "URS0000247C4D_392897-M_Cys.colored.svg",
    ]

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(self.cmd)

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestRnasep(R2dtTestCase):
    """Check that RNAse P templates work."""

    fasta_input = os.path.join("examples", "rnasep.fasta")
    test_results = os.path.join("tests", "results", "rnasep")
    precomputed_results = os.path.join("tests", "examples", "rnasep")
    cmd = f"r2dt.py rnasep draw {fasta_input} {test_results}"
    files = [
        "hits.txt",
        "URS00000A7310_29284-RNAseP_a_H_trapanicum_JB.colored.svg",
        "URS0000CBCB35_210-RNAseP_b_H_pylory_26695_JB.colored.svg",
        "URS0000EEAD19_2190-RNAseP_a_M_jannaschii_JB.colored.svg",
        "URS0001BC2932_272844-RNAseP_a_P_abyssi_JB.colored.svg",
        "URS0001BC3468_782-RNAseP_b_R_prowazekii_JB.colored.svg",
        "URS00003C82BC_186497-RNAseP_a_P_furiosus_JB.colored.svg",
        "URS00004BB8BB_511145-RNAseP_b_E_coli_JB.colored.svg",
        "URS00006A4F8D_64091-RNAseP_a_Halobacterium-NRC1_JB.colored.svg",
        "URS00006D6BE6_273075-RNAseP_a_T_acidophilum_JB.colored.svg",
        "URS00006E8172_2285-RNAseP_a_S_acidocaldarius_JB.colored.svg",
        "URS00019F4D0F_358-RNAseP_b_A_tumefaciens_JB.colored.svg",
        "URS00019F2369_1773-RNAseP_b_M_tuberculosis_JB.colored.svg",
        "URS000066E9AE_2287-RNAseP_a_S_solfataricus_JB.colored.svg",
        "URS000072E054_1095685-RNAseP_N_gonnorhoeae_JB.colored.svg",
        "URS0000637B30_1247414-RNAseP_N_gonnorhoeae_JB.colored.svg",
        "URS0000664B0C_4896-RNAseP_e_S_pombe_JB.colored.svg",
        "URS000013F331_9606-RNAseP_e_H_sapiens_3D.colored.svg",
    ]

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(self.cmd)

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestForceTemplate(R2dtTestCase):
    """Check the option that forces a sequence into a given template."""

    fasta_input = os.path.join("examples", "force")
    test_results = os.path.join("tests", "results", "force")
    precomputed_results = os.path.join("tests", "examples", "force")
    cmd = "r2dt.py draw --force_template {} {} {}"
    cases = {  # pylint: disable=duplicate-key
        "URS00000F9D45_9606": "d.5.b.E.coli",  # CRW: human 5S with E. coli 5S
        "URS0000704D22_9606": "EC_SSU_3D",  # RiboVision SSU: Human SSU with E.coli
        "URS000020CCFC_274": "EC_LSU_3D",  # RiboVision LSU: T. thermophilus with E.coli
        "URS00000A1A88_9606": "B_Thr",  # GtRNAdb: human E_Thr with B_Thr
        "URS00000A1A88_9606": "RF00005",  # GtRNAdb E_Thr using Rfam tRNA
        "URS0001BC2932_272844": "RNAseP_a_P_furiosus_JB",  # RNAse P: P. abyssi with P.furiosus
    }

    def setUp(self):
        self.delete_folder(self.test_results)
        for seq_id, model_id in self.cases.items():
            input_fasta = os.path.join(self.fasta_input, seq_id + ".fasta")
            os.system(self.cmd.format(model_id, input_fasta, self.test_results))

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        for seq_id, model_id in self.cases.items():
            filename = f"{seq_id}-{model_id}.colored.svg"
            new_file = os.path.join(self.test_results, "results", "svg", filename)
            reference_file = os.path.join(self.precomputed_results, filename)
            self.assertTrue(os.path.exists(new_file), f"File {new_file} does not exist")
            self.assertTrue(
                filecmp.cmp(new_file, reference_file),
                f"File {new_file} does not match",
            )


class TestRNAfold(R2dtTestCase):
    """Check constrained folding options."""

    fasta_input = os.path.join("examples", "constraint")
    test_results = os.path.join("tests", "results", "constraint")
    test_results_subfolder = "results/svg"
    precomputed_results = os.path.join("tests", "examples", "constraint")
    cmd = "r2dt.py draw --constraint {} {}"
    cmd2 = "r2dt.py draw --constraint --fold_type {} --force_template {} {} {}"
    fold_type_inputs = {
        "Halobacteroides_halobius1": "insertions_only",
        "Halobacteroides_halobius2": "full_molecule",
        "Halobacteroides_halobius3": "all_constraints_enforced",
    }
    output_files = {
        "Halobacteroides_halobius1-d.5.a.H.salinarum.1.colored.svg",
        "Halobacteroides_halobius2-d.5.a.H.salinarum.1.colored.svg",
        "Halobacteroides_halobius3-d.5.a.H.salinarum.1.colored.svg",
        "URS00021C62AE-RF01911.colored.svg",
        "URS0000394A9E-RF00076.colored.svg",
    }

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(
            self.cmd.format(
                os.path.join(self.fasta_input, "constraint-examples.fasta"),
                self.test_results,
            )
        )
        for seq_id, fold_type in self.fold_type_inputs.items():
            input_fasta = os.path.join(self.fasta_input, seq_id + ".fasta")
            os.system(
                self.cmd2.format(
                    fold_type, "d.5.a.H.salinarum.1", input_fasta, self.test_results
                )
            )

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestExclusions(R2dtTestCase):
    """Check constrained folding with excluded basepairs."""

    fasta_input = os.path.join(
        "examples", "constraint", "Oceanobacillus_iheyensis.fasta"
    )
    exclusion = os.path.join("examples", "constraint", "Oceanobacillus_iheyensis.txt")
    test_results = os.path.join("tests", "results", "exclusion")
    precomputed_results = os.path.join("tests", "examples", "constraint")
    cmd = f"r2dt.py draw --constraint --exclusion {exclusion} {fasta_input} {test_results}"
    output_svg = "Oceanobacillus_iheyensis-EC_SSU_3D.colored.svg"

    def setUp(self):
        self.delete_folder(self.test_results)
        os.system(self.cmd)

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        new_file = os.path.join(self.test_results, "results", "svg", self.output_svg)
        reference_file = os.path.join(self.precomputed_results, self.output_svg)
        self.assertTrue(os.path.exists(new_file), f"File {new_file} does not exist")
        self.assertTrue(
            filecmp.cmp(new_file, reference_file),
            f"File {new_file} does not match",
        )


class TestSkipRibovoreFilters(R2dtTestCase):
    """Check that the --skip_ribovore_filters option works."""

    fasta_input = os.path.join("examples", "ribovore-filters.fasta")
    test_results = os.path.join("tests", "results", "skip-ribovore-filters")
    test_results_subfolder = os.path.join("results", "svg")
    precomputed_results = os.path.join("tests", "examples", "skip-ribovore-filters")
    cmd_default = f"r2dt.py draw {fasta_input} {test_results}"
    cmd_skip = f"r2dt.py draw --skip_ribovore_filters {fasta_input} {test_results}"
    files = ["URS0000001EB3-RF00661.colored.svg"]

    def setUp(self):
        self.delete_folder(self.test_results)

    def test_default(self):
        """Check that the default command without an extra option fails."""
        os.system(self.cmd_default)
        new_file = os.path.join(
            self.test_results, self.test_results_subfolder, self.files[0]
        )
        self.assertFalse(os.path.exists(new_file), f"File {new_file} does not exist")

    def test_skip_filters(self):
        """Check that the new option works."""
        os.system(self.cmd_skip)
        self.check_examples()


if __name__ == "__main__":
    unittest.main()
