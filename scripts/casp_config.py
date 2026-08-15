"""Per-season configuration for the CASP RNA reference/model dashboard pipeline
(``casp_rank.py`` -> ``casp_fetch.py`` -> ``casp_batch.py`` -> ``casp_dashboard.py``).

Adding a new season means adding one entry to ``SEASONS`` here; the four
pipeline scripts are otherwise season-agnostic.

``target_pdb`` maps each base CASP target id (the id predictions/rankings are
keyed by) to a list of ``(state, source)`` pairs. Most targets have exactly
one experimental structure -> a single ``("", source)`` entry, and the
manifest key is just the target id. A handful of targets were solved in
multiple distinct conformations (e.g. two states of the same complex); those
list several ``(state, source)`` pairs and get one manifest entry per state
(``<target>-<state>``), all sharing the same top-N ranked model list.

``source`` is one of:
- ``rcsb(pdb_code)`` — fetch the reference from RCSB (a single structure).
- ``official(*filenames)`` — one or more files already extracted from CASP's
  own reference archive (see ``OFFICIAL_ARCHIVE_URL`` below). More than one
  filename means an *ensemble* (e.g. R1149/R1156's map-fitting uncertainty
  models): every submitted model is scored against every ensemble member and
  the best-fitting one wins — see ``casp_fetch.py``'s ``_rank_reference``.
"""

import urllib.request


def download(url, dest):
    """Idempotent download: skip if ``dest`` already exists and is non-empty."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, dest)
    return dest


def rcsb(pdb_code):
    """A reference source: a single structure fetched from RCSB."""
    return {"kind": "rcsb", "pdb": pdb_code}


def official(*filenames):
    """A reference source: one or more files from CASP's official archive."""
    return {"kind": "official", "files": list(filenames)}


# CASP16 (2024). 15-of-~44 RNA targets have a public deposited reference
# (the rest are unreleased or have no PDB deposition). No official
# CASP-hosted RNA reference archive was found for CASP16 (unlike CASP15,
# below) — RCSB is the only source. No multi-state targets in this set.
_CASP16_TARGET_PDB = {
    "R1203": [("", rcsb("8uo6"))],
    "R1205": [("", rcsb("9cfn"))],
    "R1209": [("", rcsb("9c2k"))],
    "R1211": [("", rcsb("9dcf"))],
    "R1212": [("", rcsb("9b0l"))],
    "R1242": [("", rcsb("9ely"))],
    "R1251": [("", rcsb("9mee"))],
    "R1260": [("", rcsb("9cbu"))],
    "R1261": [("", rcsb("9bzc"))],
    "R1262": [("", rcsb("9bz1"))],
    "R1263": [("", rcsb("8vqv"))],
    "R1264": [("", rcsb("8vvj"))],
    "R1285": [("", rcsb("9mcw"))],
    "R1286": [("", rcsb("9j6y"))],
    "D1273": [("", rcsb("9hio"))],
}

# CASP15 (2022). All 12 RNA targets have an official CASP reference: the
# archive below sits in the Prediction Center's public download area under
# targets/_4invitees/. Each file was verified individually — coordinates
# cross-checked bit-for-bit against the corresponding RCSB entries where
# they exist, chain naming/content confirmed by inspection — before being
# trusted here. Multi-state targets: R1136 (2 states — file names
# don't indicate which is apo/holo, so labelled v1/v2 rather than guessing),
# R1138 (young/mature RNA-origami folding intermediate, confirmed by filename),
# R1149 and R1156 (map-fitting *ensembles*, not single structures — 10 and
# 4x10 alternate models respectively, representing modelling uncertainty from
# a low-resolution cryo-EM map).
OFFICIAL_ARCHIVE_URL = (
    "https://predictioncenter.org/download_area/CASP15/targets/_4invitees/"
    "casp15_targets_RNA_4invitees.tgz"
)
# The archive also nests these two as sub-tarballs (ensembles), extracted
# alongside everything else by casp_fetch.py's archive fetch step.
_R1149_ENSEMBLE = [f"R1149_{i}.pdb" for i in range(10)]


def _r1156_conf(n):
    return [f"R1156_conformation{n}_model{i}.pdb" for i in range(1, 11)]


_CASP15_TARGET_PDB = {
    "R1107": [("", official("R1107_D_1292119758_model-annotate_P1human.pdb"))],
    "R1108": [("", official("R1108_D_1292119797_model-annotate_P1chimp.pdb"))],
    "R1116": [("", official("R1116_PV_CL_Final.pdb"))],
    "R1117": [("", official("R1117.pdb"))],
    "R1126": [("", official("R1126.pdb"))],
    "R1128": [("", official("R1128_structure_24_6wj_refined_046.pdb"))],
    "R1136": [("v1", official("R1136v1.pdb")), ("v2", official("R1136v2.pdb"))],
    "R1138": [
        ("young", official("R1138v1_structure_23v1_6hbc_young.pdb")),
        ("mature", official("R1138v2_structure_23v2_6hbc_mature.pdb")),
    ],
    "R1149": [("", official(*_R1149_ENSEMBLE))],
    "R1156": [
        ("conf1", official(*_r1156_conf(1))),
        ("conf2", official(*_r1156_conf(2))),
        ("conf3", official(*_r1156_conf(3))),
        ("conf4", official(*_r1156_conf(4))),
    ],
    "R1189": [("", official("R1189.pdb"))],
    "R1190": [("", official("R1190.pdb"))],
}


SEASONS = {
    "casp16": {
        "title": "CASP16",
        "target_pdb": _CASP16_TARGET_PDB,
        "pred_url": (
            "https://predictioncenter.org/download_area/CASP16/predictions/"
            "RNA/{target}.tar.gz"
        ),
        "rank_url": (
            "https://predictioncenter.org/download_area/CASP16/results/tables/"
            "CASP16_rna_mono.scores.csv"
        ),
        "rank_format": "casp_table",  # "Target: <id>" whitespace-column blocks
    },
    "casp15": {
        "title": "CASP15",
        "target_pdb": _CASP15_TARGET_PDB,
        "pred_url": (
            "https://predictioncenter.org/download_area/CASP15/predictions/"
            "RNA/{target}.tar.gz"
        ),
        "rank_url": "https://predictioncenter.org/download_area/CASP15/results/tables/rna.csv",
        "rank_format": "casp15_csv",  # flat CSV: target,gr_code,model,...,tm_score,...
        "official_archive_url": OFFICIAL_ARCHIVE_URL,
    },
}


def get_season(name):
    """Look up a season's config by name, or exit with the valid choices."""
    try:
        return SEASONS[name]
    except KeyError as exc:
        raise SystemExit(
            f"unknown --season {name!r}; choices: {sorted(SEASONS)}"
        ) from exc
