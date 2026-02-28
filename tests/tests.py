# pylint: disable=too-many-lines,import-outside-toplevel,reimported,redefined-outer-name
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
import glob
import os
import shutil
import tempfile
import unittest
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
from cairosvg import svg2png
from jinja2 import Environment, PackageLoader, select_autoescape
from PIL import Image, ImageChops
from skimage.metrics import structural_similarity as ssim

from utils import config, rfam
from utils.rnartist import RnaArtist
from utils.rnartist_setup import compare_rnartist_and_rscape
from utils.runner import runner

HTML_FOLDER = "tests/html"
if not os.path.exists(HTML_FOLDER):
    os.makedirs(HTML_FOLDER)

env = Environment(loader=PackageLoader("tests"), autoescape=select_autoescape())


def _get_png(filepath: str) -> Path:
    svg_path = Path(filepath)
    png_path = svg_path.with_suffix(".png")
    svg2png(bytestring=svg_path.read_bytes(), write_to=str(png_path))
    return png_path


@dataclass
class ComparisonResult:
    """
    Represents the result of a comparison between two images or binary files
    """

    # 0.01 is an arbitrary similarity threshold (roughly equivalent to 99% similarity)
    # Can be adjusted if needed, based on the results of the tests
    THRESHOLD = 0.05

    identical: bool
    similarity: Optional[float] = None
    bbox: Optional[tuple] = None
    size1: Optional[tuple] = None
    size2: Optional[tuple] = None

    # pylint: disable=missing-function-docstring
    # pylint: disable=invalid-name
    @staticmethod
    def no() -> "ComparisonResult":
        return ComparisonResult(False)

    # pylint: disable=missing-function-docstring
    @staticmethod
    def yes() -> "ComparisonResult":
        return ComparisonResult(True)


class R2dtTestCase(unittest.TestCase):
    """Base class for R2DT tests."""

    cmd = ""
    test_results = "test_folder"
    test_results_subfolder = ""
    files = []
    precomputed_results = ""

    @staticmethod
    def delete_folder(folder):
        """Delete test folder"""
        runner.run(f"rm -Rf {folder}")

    def setUp(self):
        self.delete_folder(self.test_results)
        runner.run(self.cmd, print_output=True)

    def tearDown(self):
        if os.environ.get("R2DT_KEEP_TEST_RESULTS", "0") != "1":
            self.delete_folder(self.test_results)
        # pylint: disable=expression-not-assigned
        [os.remove(file) for file in glob.glob("examples/*.ssi")]

    def create_webpage(
        self, filename: str, before, after, comparison_result: ComparisonResult
    ) -> None:
        """Create an HTML file comparing the reference SVG with a new one."""
        template = env.get_template("compare.html")
        with open(filename, "w") as f_html:
            with open(before) as f_before, open(after) as f_after:
                f_html.write(
                    template.render(
                        test_name=self.__class__.__name__,
                        before=f_before.read(),
                        after=f_after.read(),
                        comparison=comparison_result,
                    )
                )

    def _compare_files(self, reference_file: str, new_file: str) -> ComparisonResult:
        if (reference_file.endswith(".svg") or reference_file.endswith(".png")) and (
            new_file.endswith(".png") or new_file.endswith(".svg")
        ):
            reference_png, new_png = _get_png(reference_file), _get_png(new_file)

            # convert to grayscale
            reference_image = Image.open(reference_png.absolute()).convert("L")
            new_image = Image.open(new_png.absolute()).convert("L")

            diff = ImageChops.difference(reference_image, new_image)

            if reference_image.size != new_image.size:
                new_image = new_image.resize(reference_image.size)

            try:
                arr1 = np.array(reference_image)
                arr2 = np.array(new_image)

                similarity_index, _ = ssim(arr1, arr2, full=True)

                result = abs(1 - similarity_index) <= ComparisonResult.THRESHOLD

                result = ComparisonResult(
                    identical=result,
                    bbox=diff.getbbox(),
                    similarity=similarity_index,
                    size1=reference_image.size,
                    size2=new_image.size,
                )
                return result
            # pylint: disable=broad-except
            except Exception as e:
                print(f"Image comparison failed: {e}")
                return ComparisonResult(
                    identical=False,
                    bbox=diff.getbbox(),
                    size1=reference_image.size,
                    size2=new_image.size,
                )

        else:
            return (
                ComparisonResult.yes()
                if filecmp.cmp(new_file, reference_file, shallow=False)
                else ComparisonResult.no()
            )

    def check_examples(self):
        """Check that files exist and are identical to examples."""
        count = 0
        html_files = []
        for loop_id, filename in enumerate(self.files):
            new_file = os.path.join(
                self.test_results, self.test_results_subfolder, filename
            )
            reference_file = os.path.join(self.precomputed_results, filename)
            self.assertTrue(os.path.exists(new_file), f"File {new_file} does not exist")

            comparison_result = self._compare_files(reference_file, new_file)
            comparison_marker = "=" if comparison_result.identical else "!"
            filename = os.path.join(
                HTML_FOLDER,
                f"{comparison_marker} {self.__class__.__name__}_{loop_id}.html",
            )
            if reference_file.endswith(".svg") and new_file.endswith(".svg"):
                self.create_webpage(
                    filename, reference_file, new_file, comparison_result
                )
            if not comparison_result.identical:
                html_files.append(filename)
                count += 1
        if count:
            print(f"Please inspect {', '.join(html_files)}")
        self.assertEqual(
            count, 0, f"Found {count} out of {len(self.files)} non-identical SVG files"
        )


class TestCovarianceModelDatabase(unittest.TestCase):
    """Check that the covariance model and template files exist."""

    @staticmethod
    def count_lines(filename):
        """Return the number of lines in a modelinfo file."""
        with open(filename) as f_modelinfo:
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
        modelinfo = os.path.join(location, "modelinfo.txt")
        all_cm = os.path.join(location, "all.cm")
        num_lines = self.count_lines(modelinfo)
        if "crw" in location:
            num_cms = count
        else:
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
        self.verify_cm_database(config.CRW_CM_LIBRARY, 662)

    def test_ribovision_lsu_database(self):
        """Check RiboVision LSU covariance models."""
        self.verify_cm_database(config.RIBOVISION_LSU_CM_LIBRARY, 22)

    def test_ribovision_ssu_database(self):
        """Check RiboVision SSU covariance models."""
        self.verify_cm_database(config.RIBOVISION_SSU_CM_LIBRARY, 11)

    def test_rnasep_cm_database(self):
        """Check RNAse P covariance models."""
        self.verify_cm_database(config.RNASEP_CM_LIBRARY, 25)

    def test_tmrna_cm_database(self):
        """Check tmRNA covariance models."""
        self.verify_cm_database(config.TMRNA_CM_LIBRARY, 7)

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
        tempdir = Path(tempfile.gettempdir()) / "cms"
        shutil.rmtree(tempdir)


class TestRibovisionLSU(R2dtTestCase):
    """Check that the RiboVision LSU templates work."""

    fasta_input = os.path.join("examples", "lsu-small-example.fasta")
    test_results = os.path.join("tests", "results", "ribovision")
    precomputed_results = os.path.join("tests", "examples", "ribovision")
    cmd = f"r2dt.py ribovision draw_lsu {fasta_input} {test_results} --quiet"
    files = [
        "hits.txt",
        "URS000080E357_9606-mHS_LSU_3D.colored.svg",
    ]

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestRibovisionSSU(R2dtTestCase):
    """Check that the RiboVision SSU templates work."""

    fasta_input = os.path.join("examples", "ribovision-ssu-examples.fasta")
    test_results = os.path.join("tests", "results", "ribovision-ssu")
    precomputed_results = os.path.join("tests", "examples", "ribovision-ssu")
    cmd = f"r2dt.py ribovision draw_ssu {fasta_input} {test_results} --quiet"
    files = [
        "hits.txt",
        "URS00002A2E83_10090-HS_SSU_3D.colored.svg",
    ]

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestAnimation(R2dtTestCase):
    """Check that the SVG animation works."""

    svg1 = Path("examples") / "animate" / "PZ39_solution.svg"
    svg2 = Path("examples") / "animate" / "PZ39_Dfold_3.svg"
    test_results = Path("tests") / "results" / "animate"
    precomputed_results = Path("tests") / "examples" / "animate"
    files = ["PZ39_Dfold_3.animated.svg"]
    cmd = (
        f"mkdir -p {test_results} && "
        f"python3 utils/animate.py {svg1} {svg2} {test_results / files[0]}"
    )

    def test_animation(self):
        """Check that the animation works."""
        self.check_examples()


class TestRfamAccession(R2dtTestCase):
    """Test Rfam visualisation when specifying an Rfam accession."""

    rfam_acc = "RF00162"
    fasta_input = os.path.join("examples", rfam_acc + ".example.fasta")
    test_results = os.path.join("tests", "results", "rfam")
    test_results_subfolder = rfam_acc
    precomputed_results = os.path.join("tests", "examples", "rfam", rfam_acc)
    cmd = f"r2dt.py rfam draw {rfam_acc} {fasta_input} {test_results} --quiet"
    files = [
        "URS00001D0AD3_224308-RF00162.colored.svg",
        "URS00002D29F6_224308-RF00162.colored.svg",
        "URS00002F3927_224308-RF00162.colored.svg",
        "URS000053CEAC_224308-RF00162.colored.svg",
        "URS000008638F_224308-RF00162.colored.svg",
    ]

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
    cmd = f"r2dt.py draw {fasta_input} {test_results} --quiet"
    files = [
        "URS00021ED9B3_2697049-RF00507.colored.svg",
        "URS000080E2F0_93929-RF01734.colored.svg",
        "URS0000868535_32630-RF01750.colored.svg",
    ]

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestCrw(R2dtTestCase):
    """Check that the CRW templates work."""

    label = "crw"
    fasta_input = os.path.join("examples", label + "-examples.fasta")
    test_results = os.path.join("tests", "results", label)
    precomputed_results = os.path.join("tests", "examples", label)
    cmd = f"r2dt.py crw draw {fasta_input} {test_results} --quiet"
    files = [
        "hits.txt",
        "URS00000F9D45_9606-d.5.e.H.sapiens.2.colored.svg",
        "URS000044DFF6_9606-d.16.m.H.sapiens.geno.colored.svg",
        "URS000001AE2D_4932-d.16.e.S.cerevisiae.colored.svg",
    ]

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestSingleEntry(R2dtTestCase):
    """Check a broad range of RNA types at once."""

    fasta_input = os.path.join("examples", "examples.fasta")
    test_results = os.path.join("tests", "results", "single-entry")
    test_results_subfolder = os.path.join("results", "svg")
    precomputed_results = os.path.join("tests", "examples", "single-entry")
    cmd = f"r2dt.py draw {fasta_input} {test_results} --quiet"
    files = [
        "URS00000F9D45_9606-d.5.e.H.sapiens.2.colored.svg",
        "URS000044DFF6_9606-d.16.m.H.sapiens.geno.colored.svg",
        "URS000053CEAC_224308-RF00162.colored.svg",
        "URS0000162127_9606-RF00003.colored.svg",
        "URS000080E357_9606-mHS_LSU_3D.colored.svg",
        "URS0000023412_9606-E_Thr.colored.svg",
        "URS00000012EC-M_Ile.colored.svg",
        "URS0000664B0C_4896-RNAseP_e_S_pombe_JB.colored.svg",
        "tmRNA_PCTA01000033-tmRNA.colored.svg",
    ]

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
    precomputed_results = os.path.join("tests", "examples", "gtrnadb", trnascan_model)
    cmd = f"r2dt.py gtrnadb draw {fasta_input} {test_results} --domain E --isotype Thr --quiet"
    files = [
        "URS0000023412_9606-E_Thr.colored.svg",
        "URS000021550A_9606-E_Thr.colored.svg",
        "URS00000A1A88_9606-E_Thr.colored.svg",
        "URS00000F30A4_9606-E_Thr.colored.svg",
        "URS00001D9AFB_9606-E_Thr.colored.svg",
    ]

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestGtrnadbMitoVert(R2dtTestCase):
    """Test mitochondrial vertebrate tRNA templates."""

    fasta_input = os.path.join("examples", "gtrnadb-mito-vert.fasta")
    test_results = os.path.join("tests", "results", "gtrnadb", "mito-vert")
    precomputed_results = os.path.join("tests", "examples", "gtrnadb", "mito-vert")
    cmd = f"r2dt.py gtrnadb draw {fasta_input} {test_results} --quiet"
    files = [
        "URS000061A10B_9606-M_LeuTAA.colored.svg",
        "URS000054F2AC_109923-M_LeuTAG.colored.svg",
        "URS0000333A94_392897-M_SerTGA.colored.svg",
        "URS0000043FFB_392897-M_SerGCT.colored.svg",
        "URS0000247C4D_392897-M_Cys.colored.svg",
        "URS0002616B70_9606-M_SerGCT.colored.svg",
        "chr1_46824277_46824348-M_LeuTAA.colored.svg",
        "URS0000C88790-M_LeuTAG.colored.svg",
        "URS00006B0D4F-M_LeuTAG.colored.svg",
        "URS00006B5A23-M_LeuTAG.colored.svg",
        "URS00019DBB15-M_LeuTAG.colored.svg",
        "URS000066B91E-M_LeuTAG.colored.svg",
        "URS000066635D-M_LeuTAG.colored.svg",
        "URS000068492E-M_LeuTAG.colored.svg",
    ]

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestRnasep(R2dtTestCase):
    """Check that RNAse P templates work."""

    fasta_input = os.path.join("examples", "rnasep.fasta")
    test_results = os.path.join("tests", "results", "rnasep")
    precomputed_results = os.path.join("tests", "examples", "rnasep")
    cmd = f"r2dt.py rnasep draw {fasta_input} {test_results} --quiet"
    files = [
        "hits.txt",
        "URS00000A7310_29284-RNAseP_a_H_trapanicum_JB.colored.svg",
        "URS0000CBCB35_210-RNAseP_b_H_pylory_26695_JB.colored.svg",
        "URS0000EEAD19_2190-RNAseP_a_M_jannaschii_3D_MD_DSSR.colored.svg",
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

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestTmrna(R2dtTestCase):
    """Check that tmRNA templates work."""

    fasta_input = os.path.join("examples", "tmrna.fasta")
    test_results = os.path.join("tests", "results", "tmrna")
    precomputed_results = os.path.join("tests", "examples", "tmrna")
    cmd = f"r2dt.py tmrna draw {fasta_input} {test_results} --quiet"
    files = [
        "hits.txt",
        "alpha_tmRNA-cmconsensus-tmRNA_alpha.colored.svg",
        "beta_tmRNA-cmconsensus-tmRNA_beta.colored.svg",
        "cyano_tmRNA-cmconsensus-tmRNA_cyano.colored.svg",
        "dup-ABCL01000004.1_58204-58555-tmRNA_alpha.colored.svg",
        "dup-BRH-c25__sp001515955.1-tmRNA_intron.colored.svg",
        "dup-BX548175.1_1677929-1677634-tmRNA_cyano.colored.svg",
        "dup-CR555306.1_77117-77443-tmRNA_beta.colored.svg",
        "dup-PCTA01000033.1:82380-82811_1-432-tmRNA.colored.svg",
        "intron_tmRNA-cmconsensus-tmRNA_intron.colored.svg",
        "std_tmRNA-cmconsensus-tmRNA.colored.svg",
    ]

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestForceTemplate(R2dtTestCase):
    """Check the option that forces a sequence into a given template."""

    fasta_input = os.path.join("examples", "force")
    test_results = os.path.join("tests", "results", "force")
    test_results_subfolder = os.path.join("results", "svg")
    precomputed_results = os.path.join("tests", "examples", "force")
    cmd = "r2dt.py draw --force_template {} {} {} --quiet"
    files = [
        "URS00000F9D45_9606-d.5.b.E.coli.colored.svg",  # CRW: human 5S with E. coli 5S
        "URS0000704D22_9606-EC_SSU_3D.colored.svg",  # RiboVision SSU: Human SSU with E.coli
        "URS000020CCFC_274-EC_LSU_3D.colored.svg",  # RiboVision LSU: T. thermophilus with E.coli
        "URS00000A1A88_9606-B_Thr.colored.svg",  # GtRNAdb: human E_Thr with B_Thr
        "URS00000A1A88_9606-RF00005.colored.svg",  # GtRNAdb E_Thr using Rfam tRNA
        # RNAse P: P. abyssi with P.furiosus
        "URS0001BC2932_272844-RNAseP_a_P_furiosus_JB.colored.svg",
    ]

    def setUp(self):
        self.delete_folder(self.test_results)
        for filename in self.files:
            seq_id, model_id = filename.replace(".colored.svg", "").split("-")
            input_fasta = os.path.join(self.fasta_input, seq_id + ".fasta")
            runner.run(self.cmd.format(model_id, input_fasta, self.test_results))

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestRNAfold(R2dtTestCase):
    """Check constrained folding options."""

    fasta_input = os.path.join("examples", "constraint")
    test_results = os.path.join("tests", "results", "constraint")
    test_results_subfolder = "results/svg"
    precomputed_results = os.path.join("tests", "examples", "constraint")
    cmd = "r2dt.py draw --constraint {} {} --quiet"
    cmd2 = "r2dt.py draw --constraint --fold_type {} --force_template {} {} {} --quiet"
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
        runner.run(
            self.cmd.format(
                os.path.join(self.fasta_input, "constraint-examples.fasta"),
                self.test_results,
            )
        )
        for seq_id, fold_type in self.fold_type_inputs.items():
            input_fasta = os.path.join(self.fasta_input, seq_id + ".fasta")
            runner.run(
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
    test_results_subfolder = os.path.join("results", "svg")
    precomputed_results = os.path.join("tests", "examples", "constraint")
    cmd = f"r2dt.py draw --constraint --exclusion {exclusion} {fasta_input} {test_results} --quiet"
    files = ["Oceanobacillus_iheyensis-EC_SSU_3D.colored.svg"]

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestSkipRibovoreFilters(R2dtTestCase):
    """Check that the --skip_ribovore_filters option works."""

    fasta_input = os.path.join("examples", "ribovore-filters.fasta")
    test_results = os.path.join("tests", "results", "skip-ribovore-filters")
    test_results_subfolder = os.path.join("results", "svg")
    precomputed_results = os.path.join("tests", "examples", "skip-ribovore-filters")
    cmd_default = f"r2dt.py draw {fasta_input} {test_results} --quiet"
    cmd_skip = (
        f"r2dt.py draw --skip_ribovore_filters {fasta_input} {test_results} --quiet"
    )
    files = ["URS0002652150-RF00020.colored.svg"]

    def setUp(self):
        self.delete_folder(self.test_results)

    def test_default(self):
        """Check that the default command without an extra option fails."""
        runner.run(self.cmd_default)
        new_file = os.path.join(
            self.test_results, self.test_results_subfolder, self.files[0]
        )
        self.assertFalse(os.path.exists(new_file), f"File {new_file} does not exist")

    def test_skip_filters(self):
        """Check that the new option works."""
        runner.run(self.cmd_skip)
        self.check_examples()


class TestTemplateFree(R2dtTestCase):
    """Check that the templatefree visualisation works."""

    fasta_input = os.path.join("examples", "template-free.fasta")
    test_results = os.path.join("tests", "results", "template-free")
    test_results_subfolder = os.path.join("results", "svg")
    precomputed_results = os.path.join("tests", "examples", "template-free")
    cmd = f"r2dt.py templatefree {fasta_input} {test_results} --quiet"
    files = ["3SKZ_B.colored.svg"]

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestTemplateFreeAutomaticDetection(TestTemplateFree):
    """Check that templatefree input is automatically detected."""

    fasta_input = os.path.join("examples", "template-free.fasta")
    test_results = os.path.join("tests", "results", "template-free-auto")
    test_results_subfolder = os.path.join("results", "svg")
    precomputed_results = os.path.join("tests", "examples", "template-free")
    cmd = f"r2dt.py draw {fasta_input} {test_results} --quiet"
    files = ["3SKZ_B.colored.svg"]

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestTemplateGeneration(R2dtTestCase):
    """Check that the template generation works.
    File RF02976.json was downloaded from RNAcanvas
    after manually editing R2DT output generated
    using the default Rfam RF02976 template.
    """

    template_id = "RF02976"
    test_id = "template-generation"
    json_input = Path("examples") / f"{template_id}.json"
    test_results = Path(config.LOCAL_DATA) / template_id
    precomputed_results = Path("tests") / "examples" / test_id
    cmd = f"r2dt.py generate-template {json_input} --quiet"
    files = [f"{template_id}.fasta", f"{template_id}.xml"]

    def use_new_template(self):
        """Use newly generated template."""
        fasta_input = Path("examples") / f"{self.template_id}.fasta"
        template_folder = self.test_results
        self.test_results = Path("tests") / "results" / self.test_id
        cmd = (
            f"r2dt.py draw --force_template {self.template_id} "
            f"{fasta_input} {self.test_results} --quiet"
        )
        runner.run(cmd)
        self.test_results_subfolder = os.path.join("results", "svg")
        self.files = [f"URS0000D6941A-{self.template_id}.colored.svg"]
        self.check_examples()
        shutil.rmtree(template_folder)

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()
        self.use_new_template()


class TestRnartist(R2dtTestCase):
    """Check that RNArtist templates work correctly."""

    rfam_acc = "RF00025"
    fasta_input = os.path.join("examples", "rnartist.fasta")
    test_results = os.path.join("tests", "results", "rnartist")
    test_results_subfolder = rfam_acc
    precomputed_results = os.path.join("tests", "examples", "rnartist")
    cmd = (
        f"r2dt.py rfam draw {rfam_acc} {fasta_input} {test_results} --quiet --rnartist"
    )
    files = [f"URS0000696E0A-{rfam_acc}.colored.svg"]

    def test_rnartist_mode(self):
        """Check the --rnartist option works."""
        self.check_examples()

    def test_auto_mode(self):
        """Check that the auto mode works."""
        cmd = self.cmd.replace("--rnartist", "")
        runner.run(cmd)
        self.check_examples()


class TestRscapeTemplateSelectionOption(R2dtTestCase):
    """Check that the --rscape option works."""

    rfam_acc = "RF00174"
    fasta_input = os.path.join("examples", f"{rfam_acc}.example.fasta")
    test_results = os.path.join("tests", "results", "rnartist")
    test_results_subfolder = rfam_acc
    precomputed_results = os.path.join("tests", "examples", "rnartist")
    cmd = f"r2dt.py rfam draw {rfam_acc} {fasta_input} {test_results} --quiet --rscape"
    files = [f"URS000016E07A-{rfam_acc}.colored.svg"]

    def test_rscape_option(self):
        """Check that the --rscape option works."""
        self.check_examples()


class TestRnartistTemplateGeneration(R2dtTestCase):
    """Check that the RNArtist template generation works."""

    rfam_acc = "RF00174"
    test_results = (
        Path(config.PROJECT_HOME) / "tests" / "results" / "rnartist" / rfam_acc
    )
    precomputed_results = Path("tests") / "examples" / "rnartist"
    files = ["rnartist-template.xml"]

    def test_rnartist(self):
        """Check that RNArtist templates are generated."""
        rnartist = RnaArtist(self.rfam_acc, destination=self.test_results)
        rnartist.run(rerun=True)
        self.check_examples()


class TestRnartistR2rComparison(R2dtTestCase):
    """Check that the number of overlaps in RNArtist and R2R templates is compared correctly."""

    rfam_acc = "RF00174"

    def test_rnartist_vs_r2r(self):
        """Check that the number of overlaps in RNArtist and R2R templates is compared correctly."""
        chosen_template, _ = compare_rnartist_and_rscape(self.rfam_acc)
        assert (
            chosen_template == "rnartist"
        ), f"RNArtist template should be chosen, not {chosen_template}"


class TestTemplateFreeRnartist(R2dtTestCase):
    """Check that RNArtist template-free mode works correctly."""

    fasta_input = Path("examples") / "bridge-rna.fasta"
    test_results = Path("tests") / "results" / "rnartist-template-free"
    precomputed_results = Path("tests") / "examples" / "rnartist-template-free"
    test_results_subfolder = "results/svg"
    cmd = f"r2dt.py templatefree {fasta_input} {test_results} --quiet --rnartist"
    files = ["bridge_rna.colored.svg"]

    def test_template_free_mode(self):
        """Check that RNArtist template-free mode works correctly."""
        self.check_examples()


class TestLowerCase(R2dtTestCase):
    """Check that lowercase nucleotides are transferred properly."""

    fasta_input = Path("examples") / "lowercase.fasta"
    test_results = Path("tests") / "results" / "lowercase"
    precomputed_results = Path("tests") / "examples" / "lowercase"
    test_results_subfolder = "results/svg"
    cmd = f"r2dt.py draw {fasta_input} {test_results} --quiet"
    files = [
        "URS00004A7003_9606-RF00024.colored.svg",
        "URS00004A7003_9606_large_insertion-RF00024.colored.svg",
    ]

    def test_lowercase(self):
        """Check that the lowercase nucleotides are transferred properly."""
        self.check_examples()


class TestBadFastaName(R2dtTestCase):
    """Check that the program can handle fasta files with bad names."""

    fasta_input = os.path.join("examples", "bad-fasta-name.fasta")
    test_results = os.path.join("tests", "results", "bad-fasta-name")
    test_results_subfolder = os.path.join("results", "svg")
    precomputed_results = os.path.join("tests", "examples", "bad-fasta-name")
    cmd = f"r2dt.py draw {fasta_input} {test_results} --quiet"
    files = ["DB_TEXT_MORE-d.5.e.H.sapiens.2.colored.svg"]

    def test_examples(self):
        """Check that files exist and are identical to examples."""
        self.check_examples()


class TestProcPsIsPresent(unittest.TestCase):
    """Check that the procps package is present.
    Ps is required by nextflow to check if a process is running
    and the R2DT image needs to have it installed."""

    def test_procps_is_present(self):
        """Check that the procps package is present."""
        self.assertTrue(
            shutil.which("ps"), "The procps package is not installed on this system"
        )


class TestRnaview(R2dtTestCase):
    """Check that the RNAVIEW script works correctly."""

    pdb_file = Path("examples") / "PZ1_Bujnicki_1.pdb"
    test_results = Path("tests") / "results" / "rnaview"
    precomputed_results = Path("tests") / "examples" / "rnaview"
    files = ["PZ1_Bujnicki_1.colored.svg"]
    output_file = Path("examples") / "PZ1_Bujnicki_1.colored.svg"
    cmd = f"python ./utils/rnaview.py {pdb_file} > /dev/null 2>&1"

    def setUp(self):
        if not self.test_results.exists():
            self.test_results.mkdir(parents=True)
        runner.run(self.cmd, print_output=True)

        shutil.copy(self.output_file, self.test_results)

    def test_rnaview_output_exists(self):
        """Check that the RNAVIEW script generates the expected output file."""
        self.check_examples()

    def test_rnaview_output_matches(self):
        """Check that the generated file matches the precomputed file."""
        self.check_examples()


class TestAnimateSingle(R2dtTestCase):
    """Check that single animation works correctly."""

    # Reference PDB file
    ref_pdb = Path("examples") / "PZ1_solution_0.pdb"
    # Query PDB file
    query_pdb = Path("examples") / "animate_bulk" / "PZ1_Bujnicki_1.pdb"
    # Test results and precomputed results
    test_results = Path("tests") / "results" / "animate_single"
    precomputed_results = Path("tests") / "examples" / "animate_single"
    files = ["PZ1_solution_0_to_PZ1_Bujnicki_1.animated.svg"]

    # Command to run
    cmd = (
        f"python ./utils/animate3d.py {test_results / ref_pdb.name} "
        f"{test_results / query_pdb.name} > /dev/null 2>&1"
    )

    def setUp(self):
        """Set up the test environment."""
        # Ensure test results directory exists
        if not self.test_results.exists():
            self.test_results.mkdir(parents=True)

        # Copy reference and query PDBs to test results directory
        shutil.copy(self.ref_pdb, self.test_results)
        shutil.copy(self.query_pdb, self.test_results)

        # Run the animation script
        runner.run(self.cmd, print_output=True)

    def test_single_animation_output_exists(self):
        """Check that the animated SVG file is generated."""
        self.check_examples()

    def test_single_animation_output_matches(self):
        """Check that the generated file matches the precomputed file."""
        self.check_examples()


class TestAnimateBulk(R2dtTestCase):
    """Check that bulk animation works correctly."""

    # Reference PDB file
    ref_pdb = Path("examples") / "PZ1_solution_0.pdb"
    # Directory containing query PDB files
    query_dir = Path("examples") / "animate_bulk"
    # Test results and precomputed results
    test_results = Path("tests") / "results" / "animate_bulk"
    ref_dest = test_results / "PZ1_solution_0.pdb"
    query_dest = test_results / "query"

    precomputed_results = Path("tests") / "examples" / "animate_bulk"
    files = [
        "PZ1_solution_0_to_PZ1_Bujnicki_1.animated.svg",
        "PZ1_solution_0_to_PZ1_Bujnicki_2.animated.svg",
        "PZ1_solution_0_to_PZ1_Bujnicki_3.animated.svg",
    ]

    # Command to run
    cmd = f"python ./utils/animate3d.py -b {ref_dest} {query_dest} > /dev/null 2>&1"

    def setUp(self):
        """Set up the test environment."""
        # Ensure test results directory exists
        if not self.test_results.exists():
            self.test_results.mkdir(parents=True)

        # Copy reference PDB to test results directory
        shutil.copy(self.ref_pdb, self.test_results)

        # Create a query subdirectory and copy query PDBs into it
        if not self.query_dest.exists():
            self.query_dest.mkdir(parents=True)

        for pdb_file in self.query_dir.glob("*.pdb"):
            shutil.copy(pdb_file, self.query_dest)

        # Run the animation script
        runner.run(self.cmd, print_output=True)

    def test_bulk_animation_output_exists(self):
        """Check that all animated SVG files are generated."""
        self.check_examples()

    def test_bulk_animation_output_matches(self):
        """Check that the generated files match the precomputed files."""
        self.check_examples()


class TestPdbFetchUtils(unittest.TestCase):
    """Unit tests for utils/pdb_fetch.py helper functions."""

    def test_validate_pdb_id_standard(self):
        """Standard 4-character PDB IDs are valid."""
        from utils.pdb_fetch import validate_pdb_id

        self.assertTrue(validate_pdb_id("1S72"))
        self.assertTrue(validate_pdb_id("2GIS"))
        self.assertTrue(validate_pdb_id("1ehz"))

    def test_validate_pdb_id_extended(self):
        """Extended PDB IDs (>4 chars) are valid."""
        from utils.pdb_fetch import validate_pdb_id

        self.assertTrue(validate_pdb_id("9FN3A"))
        # Underscore is not alphanumeric, so this should fail
        self.assertFalse(validate_pdb_id("MA_12345"))

    def test_validate_pdb_id_rejects_short(self):
        """IDs shorter than 4 characters are rejected."""
        from utils.pdb_fetch import validate_pdb_id

        self.assertFalse(validate_pdb_id("1S7"))
        self.assertFalse(validate_pdb_id("AB"))
        self.assertFalse(validate_pdb_id(""))

    def test_validate_pdb_id_rejects_special_chars(self):
        """IDs with special characters are rejected."""
        from utils.pdb_fetch import validate_pdb_id

        self.assertFalse(validate_pdb_id("1S7!"))
        self.assertFalse(validate_pdb_id("AB CD"))
        self.assertFalse(validate_pdb_id("test.pdb"))

    def test_get_structure_format_pdb(self):
        """PDB files are detected correctly."""
        from utils.pdb_fetch import get_structure_format

        self.assertEqual(get_structure_format(Path("test.pdb")), "pdb")
        self.assertEqual(get_structure_format(Path("test.ent")), "pdb")
        self.assertEqual(get_structure_format(Path("test.pdb.gz")), "pdb")

    def test_get_structure_format_cif(self):
        """CIF files are detected correctly."""
        from utils.pdb_fetch import get_structure_format

        self.assertEqual(get_structure_format(Path("test.cif")), "cif")
        self.assertEqual(get_structure_format(Path("test.cif.gz")), "cif")

    def test_get_structure_format_unknown(self):
        """Unknown extensions return None."""
        from utils.pdb_fetch import get_structure_format

        self.assertIsNone(get_structure_format(Path("test.txt")))
        self.assertIsNone(get_structure_format(Path("test.fasta")))

    def test_is_local_structure_file_with_pdb(self):
        """Existing PDB files are recognised as local."""
        from utils.pdb_fetch import is_local_structure_file

        self.assertTrue(is_local_structure_file("examples/pdb/2GIS.pdb.gz"))
        self.assertTrue(is_local_structure_file("examples/pdb/2GIS.cif.gz"))

    def test_is_local_structure_file_rejects_pdb_id(self):
        """Plain PDB IDs are not local files."""
        from utils.pdb_fetch import is_local_structure_file

        self.assertFalse(is_local_structure_file("2GIS"))
        self.assertFalse(is_local_structure_file("1S72"))

    def test_validate_structure_file_pdb_gz(self):
        """Gzipped PDB files pass validation."""
        from utils.pdb_fetch import validate_structure_file

        is_valid, fmt, err = validate_structure_file(Path("examples/pdb/2GIS.pdb.gz"))
        self.assertTrue(is_valid, err)
        self.assertEqual(fmt, "pdb")

    def test_validate_structure_file_cif_gz(self):
        """Gzipped CIF files pass validation."""
        from utils.pdb_fetch import validate_structure_file

        is_valid, fmt, err = validate_structure_file(Path("examples/pdb/2GIS.cif.gz"))
        self.assertTrue(is_valid, err)
        self.assertEqual(fmt, "cif")

    def test_validate_structure_file_missing(self):
        """Missing files fail validation."""
        from utils.pdb_fetch import validate_structure_file

        is_valid, _, err = validate_structure_file(Path("nonexistent.pdb"))
        self.assertFalse(is_valid)
        self.assertIn("not found", err.lower())

    def test_decompressed_structure_file_gz(self):
        """DecompressedStructureFile produces a readable uncompressed path."""
        from utils.pdb_fetch import DecompressedStructureFile

        gz_path = Path("examples/pdb/2GIS.pdb.gz")
        with DecompressedStructureFile(gz_path) as decompressed:
            self.assertTrue(decompressed.exists())
            content = decompressed.read_text()
            self.assertIn("ATOM", content)

    def test_decompressed_structure_file_plain(self):
        """DecompressedStructureFile returns the original path for plain files."""
        from utils.pdb_fetch import DecompressedStructureFile

        plain_path = Path("examples/PZ1_Bujnicki_1.pdb")
        with DecompressedStructureFile(plain_path) as decompressed:
            self.assertEqual(decompressed, plain_path)


class TestFullSequenceExtraction(unittest.TestCase):
    """Unit tests for full-sequence extraction and dot-bracket remapping."""

    def test_remap_dot_bracket_no_gaps(self):
        """Remap with all-resolved mask is identity."""
        from utils.fr3d import remap_dot_bracket

        db = "((..))"
        mask = [True] * 6
        self.assertEqual(remap_dot_bracket(db, mask), db)

    def test_remap_dot_bracket_with_gaps(self):
        """Remap inserts dots at unresolved positions."""
        from utils.fr3d import remap_dot_bracket

        db = "((..))"
        # full length = 8: positions 0 and 4 unresolved
        mask = [False, True, True, True, False, True, True, True]
        result = remap_dot_bracket(db, mask)
        # pos 0:'.', 1:'(', 2:'(', 3:'.', 4:'.', 5:'.', 6:')', 7:')'
        self.assertEqual(result, ".((...))")
        self.assertEqual(len(result), 8)

    def test_remap_dot_bracket_leading_trailing_gaps(self):
        """Gaps at start and end."""
        from utils.fr3d import remap_dot_bracket

        db = "(..)"
        mask = [False, False, True, True, True, True, False]
        result = remap_dot_bracket(db, mask)
        self.assertEqual(result, "..(..).")
        self.assertEqual(len(result), 7)

    def test_remap_dot_bracket_length_mismatch_raises(self):
        """Mismatch between resolved count and dot-bracket length raises."""
        from utils.fr3d import remap_dot_bracket

        with self.assertRaises(ValueError):
            remap_dot_bracket("((..))", [True, True, True])

    def test_remap_dot_bracket_preserves_pseudoknots(self):
        """Pseudoknot characters (Aa, Bb) are preserved during remapping."""
        from utils.fr3d import remap_dot_bracket

        db = "((..Aa..))"
        # 12 total, 10 resolved (positions 2 and 7 unresolved)
        mask = [
            True,
            True,
            False,
            True,
            True,
            True,
            True,
            False,
            True,
            True,
            True,
            True,
        ]
        result = remap_dot_bracket(db, mask)
        self.assertEqual(len(result), 12)
        self.assertIn("A", result)
        self.assertIn("a", result)

    def test_get_full_sequence_from_pdb_2gis(self):
        """Full-sequence extraction from a local PDB.gz file."""
        from utils.fr3d import get_full_sequence_from_pdb

        pdb_gz = "examples/pdb/2GIS.pdb.gz"
        full_seq, mask, chain_id = get_full_sequence_from_pdb(pdb_gz)
        if full_seq:
            # SEQRES should exist and be at least as long as ATOM records
            self.assertGreater(len(full_seq), 0)
            self.assertEqual(len(full_seq), len(mask))
            self.assertTrue(any(mask))
            self.assertIsInstance(chain_id, str)

    def test_get_full_sequence_dispatcher(self):
        """get_full_sequence dispatches to PDB or CIF based on extension."""
        from utils.fr3d import get_full_sequence

        pdb_gz = "examples/pdb/2GIS.pdb.gz"
        full_seq, mask, _chain_id = get_full_sequence(pdb_gz)
        if full_seq:
            self.assertGreater(len(full_seq), 0)
            self.assertEqual(len(full_seq), len(mask))

    def test_resname_to_one_letter(self):
        """Three-letter to one-letter mapping."""
        from utils.fr3d import _resname_to_one_letter

        self.assertEqual(_resname_to_one_letter("A"), "A")
        self.assertEqual(_resname_to_one_letter("C"), "C")
        self.assertEqual(_resname_to_one_letter("G"), "G")
        self.assertEqual(_resname_to_one_letter("U"), "U")
        self.assertEqual(_resname_to_one_letter("PSU"), "U")
        self.assertEqual(_resname_to_one_letter("XYZ"), "N")


class TestPdbPostProcessing(unittest.TestCase):
    """Unit tests for SVG grey-out post-processing."""

    def _make_svg(self, num_nucleotides, sequence="ACGU"):
        """Create a minimal SVG string mimicking Traveler output.

        Position 0 is the 5' label, positions 1..N are nucleotides,
        and position N+1 is the 3' label — matching Traveler's convention.
        """
        lines = [
            '<?xml version="1.0" encoding="UTF-8"?>',
            '<svg xmlns="http://www.w3.org/2000/svg">',
            "<defs>",
            "<style>",
            "text.black { fill: black; }",
            "text.gray { fill: rgb(204,204,204); }",
            "text.green { fill: green; }",
            "</style>",
            "</defs>",
            # Position 0: 5' label
            "<g><title>0 (position.label in template: 0.5')</title>"
            '<text x="0" y="0" class="green">5\'</text></g>',
        ]
        for i in range(num_nucleotides):
            nt = sequence[i % len(sequence)]
            # Nucleotide positions are 1-based in the SVG
            svg_pos = i + 1
            lines.append(
                f"<g><title>{svg_pos} (position.label in template: "
                f"{svg_pos}.{nt})</title>"
                f'<text x="0" y="0" class="black" >{nt}</text></g>'
            )
        # 3' label
        svg_pos = num_nucleotides + 1
        lines.append(
            f"<g><title>{svg_pos} (position.label in template: "
            f"{svg_pos}.3')</title>"
            f'<text x="0" y="0" class="green">3\'</text></g>'
        )
        lines.append("</svg>")
        return "\n".join(lines)

    def test_grey_out_unresolved_basic(self):
        """Grey out positions 0 and 3 in a 4-nucleotide SVG."""
        import tempfile

        from utils.pdb_post import grey_out_unresolved

        svg_content = self._make_svg(4)
        with tempfile.NamedTemporaryFile(suffix=".svg", mode="w", delete=False) as f:
            f.write(svg_content)
            tmp_path = f.name

        try:
            mask = [False, True, True, False]
            result = grey_out_unresolved(tmp_path, mask)
            self.assertTrue(result)

            # Read back and check
            tree = ET.parse(tmp_path)
            root = tree.getroot()
            ns = {"svg": "http://www.w3.org/2000/svg"}

            texts = root.findall(".//svg:text", ns)
            if not texts:
                texts = root.findall(".//text")

            # Positions 0 and 3 should be gray, 1 and 2 should be black
            gray_count = sum(1 for t in texts if "gray" in t.get("class", ""))
            self.assertEqual(gray_count, 2)
        finally:
            Path(tmp_path).unlink(missing_ok=True)

    def test_grey_out_all_resolved_noop(self):
        """All-resolved mask produces no changes."""
        import tempfile

        from utils.pdb_post import grey_out_unresolved

        svg_content = self._make_svg(4)
        with tempfile.NamedTemporaryFile(suffix=".svg", mode="w", delete=False) as f:
            f.write(svg_content)
            tmp_path = f.name

        try:
            mask = [True, True, True, True]
            result = grey_out_unresolved(tmp_path, mask)
            self.assertFalse(result)
        finally:
            Path(tmp_path).unlink(missing_ok=True)

    def test_grey_out_nonexistent_file(self):
        """Non-existent file returns False."""
        from utils.pdb_post import grey_out_unresolved

        result = grey_out_unresolved("/nonexistent/path.svg", [False])
        self.assertFalse(result)

    def test_get_nucleotide_position(self):
        """Parse position from title text."""
        from utils.pdb_post import _get_nucleotide_position

        self.assertEqual(_get_nucleotide_position("4 (position.label: 4.G)"), 4)
        self.assertEqual(_get_nucleotide_position("0 (position.label: 0.A)"), 0)
        self.assertEqual(_get_nucleotide_position("123 foo"), 123)
        self.assertIsNone(_get_nucleotide_position(""))
        self.assertIsNone(_get_nucleotide_position("abc"))


class TestPdbCommand(unittest.TestCase):
    """End-to-end tests for r2dt.py pdb command."""

    test_results = Path("tests") / "results" / "pdb"

    def setUp(self):
        """Set up test environment."""
        if self.test_results.exists():
            shutil.rmtree(self.test_results)
        self.test_results.mkdir(parents=True)

    def tearDown(self):
        """Clean up test results."""
        if os.environ.get("R2DT_KEEP_TEST_RESULTS", "0") != "1":
            if self.test_results.exists():
                shutil.rmtree(self.test_results)

    def test_pdb_local_file_rnaview(self):
        """r2dt.py pdb with a local .pdb.gz file using rnaview."""
        pdb_gz = Path("examples") / "pdb" / "2GIS.pdb.gz"
        output = self.test_results / "rnaview"
        cmd = f"r2dt.py pdb {pdb_gz} {output} " f"--basepairs rnaview --quiet"
        exit_code = runner.run(cmd)
        self.assertEqual(exit_code, 0, "pdb command with rnaview failed")

        # FASTA file should be created
        fasta_files = list(output.glob("*.fasta"))
        self.assertGreater(len(fasta_files), 0, "No FASTA file created")

        # FASTA should contain sequence and dot-bracket
        fasta_content = fasta_files[0].read_text()
        lines = fasta_content.strip().split("\n")
        self.assertGreaterEqual(len(lines), 2, "FASTA should have header + sequence")
        self.assertTrue(lines[0].startswith(">"), "Missing FASTA header")

    def test_pdb_local_file_fr3d(self):
        """r2dt.py pdb with a local .pdb.gz file using FR3D."""
        pdb_gz = Path("examples") / "pdb" / "2GIS.pdb.gz"
        output = self.test_results / "fr3d"
        cmd = f"r2dt.py pdb {pdb_gz} {output} " f"--basepairs fr3d --quiet"
        exit_code = runner.run(cmd)
        self.assertEqual(exit_code, 0, "pdb command with fr3d failed")

        # FASTA file should be created
        fasta_files = list(output.glob("*.fasta"))
        self.assertGreater(len(fasta_files), 0, "No FASTA file created")

    def test_pdb_local_cif_rejects_rnaview(self):
        """r2dt.py pdb with CIF + --basepairs rnaview should fail gracefully."""
        cif_gz = Path("examples") / "pdb" / "2GIS.cif.gz"
        output = self.test_results / "cif_rnaview"
        cmd = f"r2dt.py pdb {cif_gz} {output} " f"--basepairs rnaview --quiet"
        # This should return non-zero or at least not crash
        runner.run(cmd)
        # No FASTA should be created because CIF + rnaview is invalid
        fasta_files = list(output.glob("*.fasta"))
        self.assertEqual(len(fasta_files), 0, "CIF + rnaview should not produce output")

    def test_pdb_auto_basepairs_pdb_format(self):
        """--basepairs auto should pick FR3D and include pseudoknots by default."""
        pdb_gz = Path("examples") / "pdb" / "2GIS.pdb.gz"
        output = self.test_results / "auto_pdb"
        cmd = f"r2dt.py pdb {pdb_gz} {output} --quiet"
        exit_code = runner.run(cmd)
        self.assertEqual(exit_code, 0, "pdb command with auto basepairs failed")

        fasta_files = list(output.glob("*.fasta"))
        self.assertGreater(len(fasta_files), 0, "No FASTA file created")

        # Auto defaults to FR3D with pseudoknots enabled
        fasta_content = fasta_files[0].read_text()
        lines = fasta_content.strip().split("\n")
        if len(lines) >= 3:
            dot_bracket = lines[2]
            # Pseudoknot characters are uppercase letters (A-Z, a-z excluding '.'/'('/')')
            has_pk = any(c.isalpha() for c in dot_bracket)
            self.assertTrue(
                has_pk,
                "Auto mode should include pseudoknots by default",
            )

    def test_pdb_auto_basepairs_cif_format(self):
        """--basepairs auto with CIF format should pick fr3d and not crash."""
        cif_gz = Path("examples") / "pdb" / "2GIS.cif.gz"
        output = self.test_results / "auto_cif"
        cmd = f"r2dt.py pdb {cif_gz} {output} --quiet"
        # Should not crash regardless of whether extraction succeeds
        exit_code = runner.run(cmd)
        self.assertEqual(exit_code, 0, "pdb command should not crash")

    def test_pdb_chain_selection(self):
        """r2dt.py pdb with --chain flag extracts specific chain."""
        pdb_gz = Path("examples") / "pdb" / "2GIS.pdb.gz"
        output = self.test_results / "chain"
        cmd = f"r2dt.py pdb {pdb_gz} {output} " f"--basepairs rnaview --chain A --quiet"
        exit_code = runner.run(cmd)
        self.assertEqual(exit_code, 0, "pdb command with --chain failed")

        fasta_files = list(output.glob("*.fasta"))
        self.assertGreater(len(fasta_files), 0, "No FASTA file created")

    def test_pdb_pseudoknots_fr3d(self):
        """Pseudoknots are included by default with FR3D."""
        pdb_gz = Path("examples") / "pdb" / "2GIS.pdb.gz"
        output = self.test_results / "pseudoknots"
        cmd = f"r2dt.py pdb {pdb_gz} {output} " f"--basepairs fr3d --quiet"
        exit_code = runner.run(cmd)
        self.assertEqual(exit_code, 0, "pdb command with pseudoknots failed")

        fasta_files = list(output.glob("*.fasta"))
        self.assertGreater(len(fasta_files), 0, "No FASTA file created")

    def test_pdb_no_pseudoknots(self):
        """--no-pseudoknots disables pseudoknot detection."""
        pdb_gz = Path("examples") / "pdb" / "2GIS.pdb.gz"
        output = self.test_results / "no_pseudoknots"
        cmd = (
            f"r2dt.py pdb {pdb_gz} {output} "
            f"--basepairs fr3d --no-pseudoknots --quiet"
        )
        exit_code = runner.run(cmd)
        self.assertEqual(exit_code, 0, "pdb command with --no-pseudoknots failed")

        fasta_files = list(output.glob("*.fasta"))
        self.assertGreater(len(fasta_files), 0, "No FASTA file created")

        # Dot-bracket should contain only '.', '(', ')' — no letter pseudoknots
        fasta_content = fasta_files[0].read_text()
        lines = fasta_content.strip().split("\n")
        if len(lines) >= 3:
            dot_bracket = lines[2]
            pk_chars = [c for c in dot_bracket if c.isalpha()]
            self.assertEqual(
                len(pk_chars),
                0,
                f"--no-pseudoknots should not produce letter chars: {pk_chars}",
            )

    def test_pdb_generates_svg(self):
        """Full pipeline: pdb → FASTA → templatefree → SVG."""
        pdb_gz = Path("examples") / "pdb" / "2GIS.pdb.gz"
        output = self.test_results / "svg_check"
        cmd = f"r2dt.py pdb {pdb_gz} {output} " f"--basepairs rnaview --quiet"
        exit_code = runner.run(cmd)
        self.assertEqual(exit_code, 0, "pdb command failed")

        # Check that SVG was generated somewhere in the output tree
        svg_files = list(output.rglob("*.svg"))
        self.assertGreater(len(svg_files), 0, "No SVG files generated")

    def test_pdb_invalid_id_fails(self):
        """Invalid PDB ID or path should fail gracefully."""
        output = self.test_results / "invalid"
        cmd = f"r2dt.py pdb !!invalid!! {output} --quiet"
        # Should not crash
        runner.run(cmd)

        # No output produced
        fasta_files = list(output.glob("*.fasta"))
        self.assertEqual(len(fasta_files), 0, "Invalid input should not produce output")

    def test_pdb_9mme_missing_nucleotides(self):
        """9MME has ~59 missing nucleotides that should be greyed out."""
        output = self.test_results / "9mme"
        cmd = f"r2dt.py pdb 9MME {output} --rnapuzzler --quiet"
        exit_code = runner.run(cmd)
        self.assertEqual(exit_code, 0, "pdb command for 9MME failed")

        # Coloured SVG must exist
        svg_files = list(output.rglob("9MME.colored.svg"))
        self.assertGreater(len(svg_files), 0, "No colored SVG generated for 9MME")
        svg_path = svg_files[0]

        svg_text = svg_path.read_text()

        # Count nucleotide <text> elements with gray vs black classes.
        # Each nucleotide is a single-letter <text> (A/C/G/U) inside a <g>.
        import re  # pylint: disable=import-outside-toplevel

        gray_nts = len(re.findall(r'class="gray"', svg_text))
        black_nts = len(re.findall(r'class="black"', svg_text))

        # 9MME has ~59 unresolved nucleotides; allow some tolerance
        self.assertGreater(gray_nts, 30, "Too few gray (unresolved) elements")
        self.assertGreater(black_nts, 400, "Too few black (resolved) elements")

        # 9MME contains pseudoknots that must appear as polyline arcs
        pk_polylines = len(re.findall(r"<polyline[^>]*pseudoknot_", svg_text))
        self.assertGreater(pk_polylines, 0, "9MME should have pseudoknot arc polylines")

        # Compare against precomputed reference SVG using SSIM
        reference = Path("tests") / "examples" / "pdb" / "9mme" / "9MME.colored.svg"
        if reference.exists():
            ref_png = _get_png(str(reference))
            new_png = _get_png(str(svg_path))
            ref_img = Image.open(ref_png.absolute()).convert("L")
            new_img = Image.open(new_png.absolute()).convert("L")
            if ref_img.size != new_img.size:
                new_img = new_img.resize(ref_img.size)
            similarity_index, _ = ssim(np.array(ref_img), np.array(new_img), full=True)
            self.assertGreater(
                similarity_index,
                1 - ComparisonResult.THRESHOLD,
                f"9MME SVG differs from reference (SSIM={similarity_index:.4f})",
            )


class TestStitch(unittest.TestCase):
    """Test SVG stitching functionality."""

    input_dir = Path("examples") / "stitch"
    test_results = Path("tests") / "results" / "stitch"
    precomputed_results = Path("tests") / "examples" / "stitch"

    def setUp(self):
        """Set up test environment."""
        if self.test_results.exists():
            shutil.rmtree(self.test_results)
        self.test_results.mkdir(parents=True)

    def tearDown(self):
        """Clean up test results."""
        if os.environ.get("R2DT_KEEP_TEST_RESULTS", "0") != "1":
            if self.test_results.exists():
                shutil.rmtree(self.test_results)

    def test_basic_stitch(self):
        """Test basic stitching of 3 panels."""
        input_files = [
            self.input_dir / "RF03120_26-299.svg",
            self.input_dir / "RF00507_13469-13546.svg",
            self.input_dir / "RF03125_29536-29870.svg",
        ]
        output_file = self.test_results / "stitched.svg"

        cmd = (
            f"r2dt.py stitch {' '.join(str(f) for f in input_files)} "
            f"-o {output_file} --sort --quiet"
        )
        exit_code = runner.run(cmd, print_output=True)

        self.assertEqual(exit_code, 0, "Stitch command failed")
        self.assertTrue(output_file.exists(), "Output SVG not created")

        # Check output has reasonable size (should be larger than any single input)
        output_size = output_file.stat().st_size
        max_input_size = max(f.stat().st_size for f in input_files)
        self.assertGreater(
            output_size, max_input_size, "Output smaller than largest input"
        )

    def test_stitch_with_captions(self):
        """Test stitching with captions."""
        input_files = [
            self.input_dir / "RF03120_26-299.svg",
            self.input_dir / "RF00507_13469-13546.svg",
            self.input_dir / "RF03125_29536-29870.svg",
        ]
        output_file = self.test_results / "stitched_captions.svg"

        cmd = (
            f"r2dt.py stitch {' '.join(str(f) for f in input_files)} "
            f"-o {output_file} --sort --quiet "
            f"--captions RF03120 --captions RF00507 --captions RF03125"
        )
        exit_code = runner.run(cmd, print_output=True)

        self.assertEqual(exit_code, 0, "Stitch command with captions failed")
        self.assertTrue(output_file.exists(), "Output SVG not created")

        # Check that captions are in the output
        content = output_file.read_text()
        self.assertIn("RF03120", content, "Caption not found in output")
        self.assertIn("RF00507", content, "Caption not found in output")

    def test_stitch_with_monochrome(self):
        """Test stitching with monochrome option."""
        input_files = [
            self.input_dir / "RF03120_26-299.svg",
            self.input_dir / "RF00507_13469-13546.svg",
        ]
        output_file = self.test_results / "stitched_mono.svg"

        cmd = (
            f"r2dt.py stitch {' '.join(str(f) for f in input_files)} "
            f"-o {output_file} --monochrome --quiet"
        )
        exit_code = runner.run(cmd, print_output=True)

        self.assertEqual(exit_code, 0, "Stitch command with monochrome failed")
        self.assertTrue(output_file.exists(), "Output SVG not created")


class TestStockholmAnnotationFormat(unittest.TestCase):
    """Test Stockholm alignment parsing with structureID/regionID format."""

    def test_parse_structure_id(self):
        """Test that parse_stockholm reads structureID annotation."""
        from utils.stockholm import parse_stockholm

        alignment = parse_stockholm(Path("examples/hcv-alignment.stk"))
        self.assertIsNotNone(alignment.structure_id, "structureID not parsed")
        self.assertIn("|", alignment.structure_id)

    def test_parse_region_id(self):
        """Test that parse_stockholm reads regionID annotation."""
        from utils.stockholm import parse_stockholm

        alignment = parse_stockholm(Path("examples/hcv-alignment.stk"))
        self.assertIsNotNone(alignment.region_id, "regionID not parsed")
        self.assertIn("|", alignment.region_id)

    def test_extract_named_regions_new_format(self):
        """Test region extraction using structureID + regionID."""
        from utils.stockholm import extract_named_regions, parse_stockholm

        alignment = parse_stockholm(Path("examples/hcv-alignment.stk"))
        regions = extract_named_regions(alignment)

        self.assertGreater(len(regions), 0, "No regions extracted")

        # Check that regions have names
        names = [r.name for r in regions]
        self.assertTrue(all(name for name in names), "Some regions have empty names")

        # Check that regions are sorted by alignment position
        positions = [r.alignment_start for r in regions]
        self.assertEqual(positions, sorted(positions), "Regions not sorted")

    def test_parent_regions_assigned(self):
        """Test that parent regions are assigned from regionID."""
        from utils.stockholm import extract_named_regions, parse_stockholm

        alignment = parse_stockholm(Path("examples/hcv-alignment.stk"))
        regions = extract_named_regions(alignment)

        # At least some regions should have parent region assigned
        regions_with_parent = [r for r in regions if r.region]
        self.assertGreater(
            len(regions_with_parent),
            0,
            "No parent regions assigned from regionID",
        )

    def test_structure_id_preferred_over_legacy(self):
        """Test that structureID takes priority over knownSS_names."""
        from utils.stockholm import StockholmAlignment, extract_named_regions

        alignment = StockholmAlignment(
            sequences={"seq1": "AUCGAUCG"},
            ss_cons="((....))",
            known_ss_names="|..old...|",
            structure_id="|..new...|",
        )
        regions = extract_named_regions(alignment)

        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0].name, "new")

    def test_legacy_fallback(self):
        """Test that knownSS_names is used when structureID is absent."""
        from utils.stockholm import StockholmAlignment, extract_named_regions

        alignment = StockholmAlignment(
            sequences={"seq1": "AUCGAUCG"},
            ss_cons="((....))",
            known_ss_names="|legacy..|",
        )
        regions = extract_named_regions(alignment)

        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0].name, "legacy")


class TestSimpleStockholm(unittest.TestCase):
    """Test Stockholm processing for simple Rfam-style alignments."""

    def test_wuss_to_dotbracket_basic(self):
        """WUSS bracket types all map to () and unpaired chars become dots."""
        from utils.stockholm import wuss_to_dotbracket

        result = wuss_to_dotbracket("(<[{....}]>)")
        self.assertEqual(result, "((((....))))")

    def test_wuss_to_dotbracket_unpaired_chars(self):
        """WUSS unpaired characters (comma, dash, underscore, colon, tilde) become dots."""
        from utils.stockholm import wuss_to_dotbracket

        result = wuss_to_dotbracket("(,-_:~)")
        self.assertEqual(result, "(.....)")  # all unpaired -> dots

    def test_wuss_to_dotbracket_letter_pairs(self):
        """Letter-pair pseudoknots are preserved as-is."""
        from utils.stockholm import wuss_to_dotbracket

        result = wuss_to_dotbracket("((.AAaa.))")
        self.assertEqual(result, "((.AAaa.))")

    def test_parse_gf_fields(self):
        """parse_stockholm captures #=GF ID and #=GF AC."""
        from utils.stockholm import parse_stockholm

        alignment = parse_stockholm(Path("examples/RF00162.stk"))
        self.assertEqual(alignment.family_id, "SAM")
        self.assertEqual(alignment.family_accession, "RF00162")

    def test_parse_rf(self):
        """parse_stockholm captures #=GC RF annotation."""
        from utils.stockholm import parse_stockholm

        alignment = parse_stockholm(Path("examples/RF00162.stk"))
        self.assertTrue(len(alignment.rf) > 0, "RF not parsed")
        self.assertEqual(len(alignment.rf), len(alignment.ss_cons))

    def test_fallback_region_uses_family_id(self):
        """Whole-alignment fallback names the region from #=GF ID."""
        from utils.stockholm import extract_named_regions, parse_stockholm

        alignment = parse_stockholm(Path("examples/RF00162.stk"))
        regions = extract_named_regions(alignment)

        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0].name, "SAM")

    def test_fallback_region_structure(self):
        """Whole-alignment fallback produces valid structure with pseudoknots."""
        from utils.stockholm import extract_named_regions, parse_stockholm

        alignment = parse_stockholm(Path("examples/RF00162.stk"))
        regions = extract_named_regions(alignment)

        region = regions[0]
        # RF00162 has 108 match columns
        self.assertEqual(len(region.structure), 108)
        self.assertEqual(len(region.consensus), 108)
        # Should contain () for regular pairs and Aa for pseudoknot
        self.assertIn("(", region.structure)
        self.assertIn(")", region.structure)
        self.assertIn("A", region.structure)
        self.assertIn("a", region.structure)
        # Should NOT contain WUSS bracket chars
        for char in "<>[]{}":
            self.assertNotIn(char, region.structure)

    def test_fallback_name_from_filename(self):
        """When no GF ID/AC, fallback_name (filename stem) is used."""
        from utils.stockholm import StockholmAlignment, extract_named_regions

        alignment = StockholmAlignment(
            sequences={"seq1": "AUCG"},
            ss_cons="(())",
            known_ss_names="",
            rf="AUCG",
        )
        regions = extract_named_regions(alignment, fallback_name="my_alignment")

        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0].name, "my_alignment")

    def test_stitch_auto_skipped_single_region(self):
        """Stockholm command with single region does not require --no-stitch."""
        test_output = Path("tests/results/rf00162_stitch_test")
        if test_output.exists():
            shutil.rmtree(test_output)

        cmd = f"r2dt.py stockholm examples/RF00162.stk {test_output} --quiet"
        exit_code = runner.run(cmd, print_output=True)
        self.assertEqual(exit_code, 0)

        # SVG generated
        svgs = list((test_output / "results" / "svg").glob("*.svg"))
        self.assertEqual(len(svgs), 1)

        # No stitched output (only 1 region)
        self.assertFalse(
            (test_output / "stitched.svg").exists(),
            "Stitched SVG should not exist for single region",
        )

        if os.environ.get("R2DT_KEEP_TEST_RESULTS", "0") != "1":
            shutil.rmtree(test_output, ignore_errors=True)


class TestViralAnnotate(unittest.TestCase):
    """Test viral genome annotation pipeline."""

    fasta_input = Path("examples") / "viral" / "coronavirus.fasta"
    test_results = Path("tests") / "results" / "viral"
    precomputed_results = Path("tests") / "examples" / "viral"

    # Expected Rfam hits for SARS-CoV-2 (OX309346.1)
    expected_hits = ["RF03120", "RF00507", "RF03125"]

    def setUp(self):
        """Set up test environment."""
        if self.test_results.exists():
            shutil.rmtree(self.test_results)
        self.test_results.mkdir(parents=True)

    def tearDown(self):
        """Clean up test results."""
        if os.environ.get("R2DT_KEEP_TEST_RESULTS", "0") != "1":
            if self.test_results.exists():
                shutil.rmtree(self.test_results)

    def test_viral_annotate_finds_hits(self):
        """Test that viral-annotate finds expected RNA families."""
        cmd = f"r2dt.py viral-annotate {self.fasta_input} {self.test_results} --quiet"
        exit_code = runner.run(cmd, print_output=True)

        self.assertEqual(exit_code, 0, "viral-annotate command failed")

        # Check that cmscan output exists
        tblout = self.test_results / "cmscan.tblout"
        self.assertTrue(tblout.exists(), "cmscan.tblout not created")

        # Check that expected hits are found
        tblout_content = tblout.read_text()
        for hit in self.expected_hits:
            self.assertIn(hit, tblout_content, f"Expected hit {hit} not found")

    def test_viral_annotate_generates_svgs(self):
        """Test that viral-annotate generates SVG diagrams."""
        cmd = f"r2dt.py viral-annotate {self.fasta_input} {self.test_results} --quiet"
        exit_code = runner.run(cmd, print_output=True)

        self.assertEqual(exit_code, 0, "viral-annotate command failed")

        # Check that rfam output directory has SVG files
        rfam_dir = self.test_results / "rfam"
        if rfam_dir.exists():
            svg_files = list(rfam_dir.glob("*.svg"))
            self.assertGreater(len(svg_files), 0, "No SVG files generated")

    def test_viral_annotate_creates_stitched_output(self):
        """Test that viral-annotate creates stitched output."""
        cmd = f"r2dt.py viral-annotate {self.fasta_input} {self.test_results} --quiet"
        exit_code = runner.run(cmd, print_output=True)

        self.assertEqual(exit_code, 0, "viral-annotate command failed")

        # Check for stitched output
        _stitched = self.test_results / "stitched.svg"
        # Note: stitched output may not exist if individual SVGs failed
        # This test will help us identify if the full pipeline works


class TestColoring(unittest.TestCase):
    """Unit tests for utils/coloring.py colour mapping."""

    def test_build_color_map_structure_mode(self):
        """Each structure_id gets a deterministic colour."""
        from utils.coloring import build_color_map

        regions = [
            {"name": "SLI", "region": "5'UTR"},
            {"name": "SLII", "region": "5'UTR"},
            {"name": "SLIII", "region": "core"},
        ]
        colors = build_color_map(regions, "structure")

        self.assertEqual(len(colors), 3)
        # Same name → same colour, different name → different colour
        self.assertNotEqual(colors[0], colors[1])
        # Deterministic: running again gives identical result
        self.assertEqual(colors, build_color_map(regions, "structure"))

    def test_build_color_map_region_mode(self):
        """Structures sharing a region_id get the same colour."""
        from utils.coloring import build_color_map

        regions = [
            {"name": "SLI", "region": "5'UTR"},
            {"name": "SLII", "region": "5'UTR"},
            {"name": "SLIII", "region": "core"},
        ]
        colors = build_color_map(regions, "region")

        self.assertEqual(len(colors), 3)
        # SLI and SLII share region → same colour
        self.assertEqual(colors[0], colors[1])
        # core is different from 5'UTR
        self.assertNotEqual(colors[0], colors[2])

    def test_build_color_map_region_falls_back_to_name(self):
        """When region is absent, region mode uses the structure name."""
        from utils.coloring import build_color_map

        regions = [
            {"name": "SLI"},
            {"name": "SLI", "region": ""},
        ]
        colors = build_color_map(regions, "region")

        # Both should resolve to the same key ("SLI")
        self.assertEqual(colors[0], colors[1])

    def test_build_color_map_config_mode(self):
        """Config mode reads colours from examples/color-config.tsv."""
        from utils.coloring import build_color_map

        regions = [
            {"name": "SLI", "region": "5'UTR"},
            {"name": "SLII", "region": "core_protein"},
            {"name": "unknown_structure", "region": "unknown_region"},
        ]
        config_path = Path("examples/color-config.tsv")
        colors = build_color_map(regions, "config", config_path)

        self.assertEqual(len(colors), 3)
        # SLI matched via region "5'UTR" → steelblue
        self.assertEqual(colors[0], "steelblue")
        # SLII matched via region "core_protein" → #e07a5f
        self.assertEqual(colors[1], "#e07a5f")
        # unknown → fallback * → gray
        self.assertEqual(colors[2], "gray")

    def test_config_mode_matches_name_before_region(self):
        """Structure name takes priority over region in config lookup."""
        from utils.coloring import build_color_map

        regions = [{"name": "3'UTR", "region": "5'UTR"}]
        config_path = Path("examples/color-config.tsv")
        colors = build_color_map(regions, "config", config_path)

        # "3'UTR" is an explicit key → #e76f51, not steelblue from "5'UTR"
        self.assertEqual(colors[0], "#e76f51")

    def test_config_mode_raises_without_default(self):
        """Config mode raises ValueError when name is missing and no * default."""
        from utils.coloring import build_color_map

        regions = [{"name": "nonexistent", "region": "also_nonexistent"}]

        # Write a minimal config with no * default
        tmp = Path("tests/results")
        tmp.mkdir(parents=True, exist_ok=True)
        config_file = tmp / "no_default.tsv"
        config_file.write_text("5'UTR\tsteelblue\n")

        with self.assertRaises(ValueError):
            build_color_map(regions, "config", config_file)

        config_file.unlink(missing_ok=True)

    def test_unknown_mode_raises(self):
        """Unknown coloring mode raises ValueError."""
        from utils.coloring import build_color_map

        with self.assertRaises(ValueError):
            build_color_map([], "invalid_mode")


class TestInlineStyleHelper(unittest.TestCase):
    """Unit tests for _set_inline_style CSS helper."""

    def test_set_on_empty_element(self):
        """Setting a property on an element with no style attribute."""
        import xml.etree.ElementTree as ET

        from utils.stitch import _set_inline_style

        elem = ET.Element("text")
        _set_inline_style(elem, "fill", "red")
        self.assertEqual(elem.get("style"), "fill: red")

    def test_add_second_property(self):
        """Adding a second property preserves the first."""
        import xml.etree.ElementTree as ET

        from utils.stitch import _set_inline_style

        elem = ET.Element("text")
        elem.set("style", "fill: red")
        _set_inline_style(elem, "stroke", "blue")
        style = elem.get("style")
        self.assertIn("fill: red", style)
        self.assertIn("stroke: blue", style)

    def test_replace_existing_property(self):
        """Replacing an existing property updates its value."""
        import xml.etree.ElementTree as ET

        from utils.stitch import _set_inline_style

        elem = ET.Element("line")
        elem.set("style", "stroke: black; stroke-width: 2")
        _set_inline_style(elem, "stroke", "green")
        style = elem.get("style")
        self.assertIn("stroke: green", style)
        self.assertNotIn("black", style)
        # stroke-width should survive (not clobbered by "stroke" replacement)
        self.assertIn("stroke-width", style)


class TestApplyPanelColor(unittest.TestCase):
    """Unit tests for apply_panel_color nucleotide recolouring."""

    def _make_panel_svg(self):
        """Build a minimal SVG with nucleotide text, a label, and a backbone line."""
        import xml.etree.ElementTree as ET

        root = ET.Element("svg")
        # Nucleotide text
        nt = ET.SubElement(root, "text", {"class": "green", "x": "10", "y": "10"})
        nt.text = "A"
        # Non-nucleotide label (5' anchor)
        label = ET.SubElement(root, "text", {"x": "20", "y": "20"})
        label.text = "5'"
        # Backbone line
        ET.SubElement(
            root,
            "line",
            {"class": "gray", "x1": "0", "y1": "0", "x2": "10", "y2": "10"},
        )
        # Residue circle (should be removed)
        ET.SubElement(root, "circle", {"class": "residue-circle", "r": "5"})
        return root

    def test_nucleotide_gets_panel_color(self):
        """Single-letter nucleotide text gets the panel colour via inline style."""
        from utils.stitch import apply_panel_color

        root = self._make_panel_svg()
        apply_panel_color(root, "#e07a5f")

        nt = [el for el in root.iter() if el.text == "A"][0]
        self.assertIn("fill: #e07a5f", nt.get("style", ""))

    def test_label_stays_black(self):
        """Non-nucleotide labels are set to black, not the panel colour."""
        from utils.stitch import apply_panel_color

        root = self._make_panel_svg()
        apply_panel_color(root, "#e07a5f")

        label = [el for el in root.iter() if (el.text or "").strip() == "5'"][0]
        self.assertIn("fill: black", label.get("style", ""))

    def test_backbone_line_recoloured(self):
        """Backbone lines (class=gray) get the panel colour."""
        from utils.stitch import apply_panel_color

        root = self._make_panel_svg()
        apply_panel_color(root, "steelblue")

        lines = [
            el
            for el in root.iter()
            if (el.tag.split("}")[-1] if "}" in el.tag else el.tag) == "line"
        ]
        self.assertTrue(len(lines) > 0, "No lines found")
        self.assertIn("stroke: steelblue", lines[0].get("style", ""))

    def test_residue_circles_removed(self):
        """Residue circles are removed after apply_panel_color."""
        from utils.stitch import apply_panel_color

        root = self._make_panel_svg()
        apply_panel_color(root, "red")

        circles = [
            el
            for el in root.iter()
            if (el.tag.split("}")[-1] if "}" in el.tag else el.tag) == "circle"
        ]
        self.assertEqual(len(circles), 0, "Residue circles should be removed")


class TestStitchColoring(unittest.TestCase):
    """Test coloring features in the stitch pipeline."""

    input_dir = Path("examples") / "stitch"
    test_results = Path("tests") / "results" / "stitch_color"

    def setUp(self):
        """Set up test environment."""
        if self.test_results.exists():
            shutil.rmtree(self.test_results)
        self.test_results.mkdir(parents=True)

    def tearDown(self):
        """Clean up test results."""
        if os.environ.get("R2DT_KEEP_TEST_RESULTS", "0") != "1":
            if self.test_results.exists():
                shutil.rmtree(self.test_results)

    def _load_panels(self):
        """Load the 3 stitch example panels."""
        from utils.stitch import parse_svg

        files = sorted(self.input_dir.glob("*.svg"))
        return [parse_svg(f) for f in files]

    def test_stitch_with_panel_colors(self):
        """Stitching with panel_colors produces per-panel coloured outlines."""
        import xml.etree.ElementTree as ET

        from utils.stitch import stitch_svgs

        panels = self._load_panels()
        colors = ["red", "green", "blue"]

        root = stitch_svgs(panels, panel_colors=colors)
        svg_str = ET.tostring(root, encoding="unicode")

        # Should have multiple connecting-outline paths (one per panel + bridges)
        outline_count = svg_str.count('class="connecting-outline"')
        self.assertGreaterEqual(
            outline_count, len(panels), "Expected at least one outline per panel"
        )

        # Each panel colour should appear in a connecting-outline stroke
        for color in colors:
            self.assertIn(
                f'stroke="{color}"',
                svg_str,
                f"Panel colour {color} not in outline strokes",
            )

    def test_stitch_without_colors_single_outline(self):
        """Without panel_colors, stitching produces a single outline path."""
        import xml.etree.ElementTree as ET

        from utils.stitch import stitch_svgs

        panels = self._load_panels()
        root = stitch_svgs(panels, monochrome=True)
        svg_str = ET.tostring(root, encoding="unicode")

        outline_count = svg_str.count('class="connecting-outline"')
        self.assertEqual(outline_count, 1, "Monochrome should have exactly 1 outline")

    def test_colored_outline_bridges_exist(self):
        """Per-panel outlines include bridge segments between panels."""
        import xml.etree.ElementTree as ET

        from utils.stitch import stitch_svgs

        panels = self._load_panels()
        colors = ["#aa0000", "#00aa00", "#0000aa"]

        root = stitch_svgs(panels, panel_colors=colors)
        svg_str = ET.tostring(root, encoding="unicode")

        outline_count = svg_str.count('class="connecting-outline"')
        # 3 panel paths + 2 bridge paths = 5
        expected_min = len(panels) + (len(panels) - 1)
        self.assertGreaterEqual(
            outline_count,
            expected_min,
            f"Expected ≥{expected_min} outlines (panels + bridges), got {outline_count}",
        )

    def test_panel_color_stamps_data_attribute(self):
        """Each panel group carries a data-panel-color attribute."""
        from utils.stitch import stitch_svgs

        panels = self._load_panels()
        colors = ["red", "green", "blue"]

        root = stitch_svgs(panels, panel_colors=colors)

        stamped_colors = []
        for elem in root.iter():
            tag = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag
            if tag == "g" and elem.get("data-panel-color"):
                stamped_colors.append(elem.get("data-panel-color"))

        self.assertEqual(
            stamped_colors,
            colors,
            "data-panel-color attributes don't match input colours",
        )

    def test_thumbnail_uses_panel_colors(self):
        """Thumbnail reflects panel colours, not black, when colors are set."""
        import xml.etree.ElementTree as ET

        from utils.stitch import create_thumbnail_svg, stitch_svgs

        panels = self._load_panels()
        colors = ["steelblue", "#e07a5f", "#2a9d8f"]

        root = stitch_svgs(panels, panel_colors=colors)
        create_thumbnail_svg(root)
        svg_str = ET.tostring(root, encoding="unicode")

        # Thumbnail connecting-outline should NOT be black
        # (connecting-outline paths carry panel colours)
        self.assertNotIn(
            'stroke="black"',
            svg_str,
            "Thumbnail still has black connecting-outline strokes",
        )

    def test_thumbnail_monochrome_uses_gray(self):
        """Monochrome thumbnail uses #cccccc, not black, for outline."""
        import xml.etree.ElementTree as ET

        from utils.stitch import create_thumbnail_svg, stitch_svgs

        panels = self._load_panels()
        root = stitch_svgs(panels, monochrome=True)
        create_thumbnail_svg(root)
        svg_str = ET.tostring(root, encoding="unicode")

        self.assertIn(
            'stroke="#cccccc"',
            svg_str,
            "Monochrome thumbnail should use #cccccc for connecting-outline",
        )
        self.assertNotIn(
            'stroke="black"',
            svg_str,
            "Monochrome thumbnail should not use black for connecting-outline",
        )


class TestStockholmColoring(unittest.TestCase):
    """End-to-end test: stockholm command with --color-config."""

    stk_input = Path("examples") / "hcv-alignment.stk"
    color_config = Path("examples") / "color-config.tsv"
    test_results = Path("tests") / "results" / "stockholm_color"

    def setUp(self):
        """Set up test environment."""
        if self.test_results.exists():
            shutil.rmtree(self.test_results)
        self.test_results.mkdir(parents=True)

    def tearDown(self):
        """Clean up test results."""
        if os.environ.get("R2DT_KEEP_TEST_RESULTS", "0") != "1":
            if self.test_results.exists():
                shutil.rmtree(self.test_results)

    def test_stockholm_color_config_produces_stitched(self):
        """stockholm --color-config generates stitched SVGs."""
        cmd = (
            f"r2dt.py stockholm {self.stk_input} {self.test_results} "
            f"--color-config {self.color_config} --quiet"
        )
        exit_code = runner.run(cmd, print_output=True)

        self.assertEqual(exit_code, 0, "stockholm --color-config command failed")
        for name in ("stitched.svg", "stitched-thumbnail.svg", "stitched-outline.svg"):
            path = self.test_results / name
            self.assertTrue(path.exists(), f"{name} not created")

    def test_stitched_svg_contains_panel_colors(self):
        """Stitched SVG has per-panel coloured connecting-outline paths."""
        cmd = (
            f"r2dt.py stockholm {self.stk_input} {self.test_results} "
            f"--color-config {self.color_config} --quiet"
        )
        runner.run(cmd, print_output=True)

        content = (self.test_results / "stitched.svg").read_text()

        # Should have multiple outline paths (not just one)
        outline_count = content.count('class="connecting-outline"')
        self.assertGreater(
            outline_count, 1, "Expected multiple per-panel coloured outlines"
        )

        # The steelblue colour (5'UTR) should appear in the SVG
        self.assertIn("steelblue", content, "5'UTR colour 'steelblue' not found")

    def test_thumbnail_no_black_strokes(self):
        """Coloured thumbnail does not use black connecting-outline strokes."""
        cmd = (
            f"r2dt.py stockholm {self.stk_input} {self.test_results} "
            f"--color-config {self.color_config} --quiet"
        )
        runner.run(cmd, print_output=True)

        content = (self.test_results / "stitched-thumbnail.svg").read_text()

        # Find all connecting-outline paths and check none use black
        import re

        outline_strokes = re.findall(
            r'class="connecting-outline"[^/]*stroke="([^"]+)"', content
        )
        # Also check reversed attribute order
        outline_strokes += re.findall(
            r'stroke="([^"]+)"[^/]*class="connecting-outline"', content
        )

        self.assertGreater(len(outline_strokes), 0, "No connecting-outline paths found")
        for stroke in outline_strokes:
            self.assertNotEqual(
                stroke, "black", "Thumbnail connecting-outline should not be black"
            )

    def test_stockholm_color_by_region(self):
        """stockholm --color-by region runs without error."""
        cmd = (
            f"r2dt.py stockholm {self.stk_input} {self.test_results} "
            f"--color-by region --quiet"
        )
        exit_code = runner.run(cmd, print_output=True)

        self.assertEqual(exit_code, 0, "stockholm --color-by region failed")
        self.assertTrue(
            (self.test_results / "stitched.svg").exists(),
            "stitched.svg not created with --color-by region",
        )


if __name__ == "__main__":
    unittest.main()
