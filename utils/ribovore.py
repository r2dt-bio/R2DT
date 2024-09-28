"""
A Python class that parses the output of RiboVore and decides which hits to keep.
"""

from dataclasses import dataclass
from pathlib import Path

from utils.rfam_threshold import RfamThresholds

MIN_GA = 20  # minimum acceptable score for Rfam models


# pylint: disable=too-many-instance-attributes, line-too-long
@dataclass
class RibovoreHit:
    """A dataclass to store information about a hit from RiboVore output."""

    idx: int  # index of sequence in input sequence file
    target: str  # name of target sequence
    pass_or_fail: str  # PASS or FAIL
    length: int  # length of target sequence (nt)
    num_families: int  # number of different families detected in sequence
    family: str  # name of family the best-scoring model to this sequence belongs to
    domain: str  # name of domain the best-scoring model to this sequence belongs to
    model: str  # name of best-scoring model
    strand: str  # strand ('plus' or 'minus') of best-scoring hit
    num_hits: int  # number of hits to best model on strand of best hit (no score threshold enforced)
    tscore: float  # summed bit scores of all hits between best model and this sequence (no score threshold enforced)
    bscore: float  # bit score of best-scoring hit between best model and this sequence (above threshold)
    s_nt: float  # summed bit scores of all hits divided by length of the sequence
    best_evalue: float  # E-value of best-scoring hit to this sequence
    tcov: float  # fraction of target sequence included in all (non-overlapping) hits to the best-scoring model
    bcov: float  # fraction of target sequence included in single best-scoring hit
    bfrom: int  # start position in the sequence of best-scoring hit
    bto: int  # stop position in the sequence of best-scoring hit
    mfrom: int  # start position in the model of best-scoring hit
    mto: int  # stop position in the model of best-scoring hit
    scdiff: float  # difference in score from classification stage between summed score of hits to best model and summed scores of hits to second best model
    scd_nt: float  # score difference per position: 'scdiff' value divided by total length of all hits to best model
    second_model: str  # name of second best-scoring model
    second_tscore: float  # summed bit scores of all hits between second-best model and this sequence (no score threshold enforced)
    unexpected_features: str  # unexpected/unusual features of sequence

    def __str__(self):
        return f"{self.target} {self.pass_or_fail}"


class Ribovore:
    """A class to parse the output of RiboVore and decide which hits to keep."""

    def __init__(self, ribovore_output: Path, skip_ribovore_filters: bool = False):
        self.ribovore_output = ribovore_output
        self.skip_ribovore_filters = skip_ribovore_filters
        self.all_hits = self.parse_ribovore_output()
        self.hits = self.filter_hits()

    # pylint: disable=too-many-locals
    def parse_ribovore_output(self):
        """
        #idx  target                  p/f  length  #fm  fam  domain    model                          strnd  #ht  tscore  bscore  s/nt   bevalue   tcov   bcov   bfrom     bto  mfrom    mto  scdiff  scd/nt  model                          tscore  unexpected_features
        #---  ---------------------  ----  ------  ---  ---  --------  -----------------------------  -----  ---  ------  ------  ----  --------  -----  -----  ------  ------  -----  -----  ------  ------  -----------------------------  ------  -------------------
        1     URS0000005270_10090    FAIL     157    0  -    -         -                              -        -       -       -     -         -      -      -       -       -      -      -       -       -  -                                   -  *NoHits;
        """
        ribovore_hits = []
        with open(self.ribovore_output) as f_ribovore:
            for line in f_ribovore:
                if line.startswith("#"):
                    continue
                parts = line.strip().split()
                idx = int(parts[0])
                target = parts[1]
                pass_or_fail = parts[2]
                length = int(parts[3])
                num_families = int(parts[4])
                family = parts[5]
                domain = parts[6]
                model = parts[7]
                strand = parts[8]
                num_hits = 0 if parts[9] == "-" else int(parts[9])
                tscore = 0 if parts[10] == "-" else float(parts[10])
                bscore = 0 if parts[11] == "-" else float(parts[11])
                s_nt = 0 if parts[12] == "-" else float(parts[12])
                best_evalue = 0 if parts[13] == "-" else float(parts[13])
                tcov = 0 if parts[14] == "-" else float(parts[14])
                bcov = 0 if parts[15] == "-" else float(parts[15])
                bfrom = 0 if parts[16] == "-" else int(parts[16])
                bto = 0 if parts[17] == "-" else int(parts[17])
                mfrom = 0 if parts[18] == "-" else int(parts[18])
                mto = 0 if parts[19] == "-" else int(parts[19])
                scdiff = 0 if parts[20] == "-" else float(parts[20])
                scd_nt = 0 if parts[21] == "-" else float(parts[21])
                second_model = None if parts[22] == "-" else parts[22]
                second_tscore = 0 if parts[23] == "-" else float(parts[23])
                unexpected_features = parts[-1]

                ribovore_hit = RibovoreHit(
                    idx,
                    target,
                    pass_or_fail,
                    length,
                    num_families,
                    family,
                    domain,
                    model,
                    strand,
                    num_hits,
                    tscore,
                    bscore,
                    s_nt,
                    best_evalue,
                    tcov,
                    bcov,
                    bfrom,
                    bto,
                    mfrom,
                    mto,
                    scdiff,
                    scd_nt,
                    second_model,
                    second_tscore,
                    unexpected_features,
                )

                ribovore_hits.append(ribovore_hit)
        return ribovore_hits

    def filter_hits(self):
        """
        Filter the hits based on the Rfam model-specific gathering thresholds.
        """
        rfam_thresholds = RfamThresholds()
        filtered_hits = []
        for hit in self.all_hits:
            if "NoHits" in hit.unexpected_features:
                continue
            if self.skip_ribovore_filters:
                filtered_hits.append(hit)
            else:
                if "NoHits" in hit.unexpected_features:
                    continue
                if "MinusStrand" in hit.unexpected_features:
                    continue
                if hit.tcov < 0.2:
                    continue
                # For Rfam models, check if the hit passes the gathering threshold
                # Note that Ribovore scores are HMM bit scores, not CM scores like Rfam GA
                rfam_ga = rfam_thresholds.get_threshold(hit.model)
                if rfam_ga is None:  # not an Rfam model, keep the hit
                    filtered_hits.append(hit)
                else:  # Rfam model, check if hit passes the gathering threshold
                    if hit.tscore >= MIN_GA:
                        filtered_hits.append(hit)
        return filtered_hits
