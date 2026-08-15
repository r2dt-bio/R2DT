"""Compare-viewer emission for the ``pdb`` command: reference vs model
2D+3D reports, plus the shared viewer-asset copying.

Split out of r2dt.py mechanically; the CLI commands remain there.
"""

# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
# pylint: disable=too-many-branches,too-many-statements,duplicate-code

import json
import re
import shutil
from pathlib import Path

import click  # pylint: disable=import-error
from rich import print as rprint

from utils import lbn_export, viewer_export, viewer_html

# Set by r2dt.py after it defines the ``templatefree`` click command: these
# helpers re-run the 2D layout via ``ctx.invoke(templatefree, ...)``, and the
# command object lives in r2dt.py, which imports this module. Injecting the
# attribute avoids the circular import.
templatefree = None  # pylint: disable=invalid-name


def ensure_cif(path, workdir, quiet):
    """Return an mmCIF path for ``path``, converting a PDB via BioPython.

    ``utils.multichain`` reads structures through FR3D's mmCIF reader, so a
    predicted model supplied as ``.pdb`` is converted first.
    """
    if path.suffix.lower() == ".cif":
        return path
    # pylint: disable=import-outside-toplevel
    import warnings

    from Bio.PDB import MMCIFIO, PDBParser

    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / f"{path.stem}.cif"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = PDBParser(QUIET=True).get_structure(path.stem, str(path))
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(str(out))
    if not quiet:
        rprint(f"[dim]Converted model {path.name} -> {out.name}[/dim]")
    return out


# CASP16 assessment colour scheme (Das et al.): the experimental structure is
# coloured by molecule type — RNA green, DNA orange, ligand red, protein yellow —
# and the best predicted model is overlaid in dark blue.
_CASP_RNA_GREEN = {"r": 0, "g": 166, "b": 81}
_CASP_DNA_ORANGE = {"r": 247, "g": 148, "b": 30}
_CASP_LIGAND_RED = {"r": 237, "g": 28, "b": 36}
_CASP_PROTEIN_YELLOW = {"r": 255, "g": 204, "b": 0}
_CASP_MODEL_BLUE = {"r": 26, "g": 58, "b": 140}


def layout_multichain_structure(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    ctx,
    *,
    structure_id: str,
    sequence: str,
    dot_bracket: str,
    boundaries,
    nested_pairs,
    crossing_pairs,
    layout_dir: Path,
    rnapuzzler_flag: bool,
    quiet: bool,
):
    """Run templatefree on a multi-chain FASTA; return (colored_json, colored_svg) or None."""
    # pylint: disable=import-outside-toplevel
    from utils import multichain

    layout_dir = Path(layout_dir)
    layout_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = layout_dir / f"{structure_id}.fasta"
    fasta_path.write_text(
        f">{structure_id}\n{sequence}\n{dot_bracket}\n", encoding="utf-8"
    )
    results_folder = layout_dir / "results"
    ctx.invoke(
        templatefree,
        fasta_input=str(fasta_path),
        output_folder=str(results_folder),
        rnartist=False,
        rscape=False,
        rnapuzzler_flag=rnapuzzler_flag,
        quiet=quiet,
    )
    svg_dir = results_folder / "results" / "svg"
    json_dir = results_folder / "results" / "json"
    svg_path = next(
        (
            p
            for p in (
                svg_dir / f"{structure_id}.colored.svg",
                svg_dir / f"{structure_id}.svg",
            )
            if p.exists()
        ),
        None,
    )
    json_path = json_dir / f"{structure_id}.colored.json"
    if svg_path is None or not json_path.exists():
        return None
    multichain.postprocess_combined_svg(
        str(svg_path),
        boundaries,
        nested_pairs,
        crossing_pairs,
        quiet=quiet,
    )
    return json_path, svg_path


def emit_compare_viewer(  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-statements,too-many-branches
    ctx,
    source_structure_path,
    actual_format,
    structure_id,
    output_path,
    result,
    simulate_seed,
    quiet,
    model_file=None,
    model_chains=None,
    score_chains=None,
    chain_views=None,
    rnapuzzler_flag=False,
):
    """Assemble the interactive 3-panel compare page.

    Reference 2D and the default model 2D panel share the reference's combined
    layout (approach B). A second model layout (``model-own/``) is also
    generated from the model's own structure so the viewer can switch on
    demand. Each panel's data lives under ``ref/`` / ``model/`` / ``model-own/``
    next to ``index.html``.
    """
    # pylint: disable=import-outside-toplevel
    from utils import multichain

    results_folder = output_path / "results"
    colored_json = results_folder / "results" / "json" / f"{structure_id}.colored.json"
    colored_svg = results_folder / "results" / "svg" / f"{structure_id}.colored.svg"
    if not colored_json.exists():
        if not quiet:
            rprint("[yellow]No colored JSON; cannot build compare viewer.[/yellow]")
        return

    n = len(result.sequence)

    # score_offsets[k] = the absolute position in result's (possibly widened)
    # label space that the model's k-th residue (in its own, narrower, local
    # order) corresponds to. Identity when unwidened, so every computation
    # below that uses it is a no-op in the common case.
    if score_chains:
        by_chain = {cid: (start, end) for cid, start, end in result.boundaries}
        score_offsets = [pos for cid in score_chains for pos in range(*by_chain[cid])]
    else:
        score_offsets = list(range(n))
    score_positions = set(score_offsets)
    score_n = len(score_offsets)

    def _remap_pairs(pairs):
        return [(score_offsets[i], score_offsets[j]) for i, j in pairs]

    def _remap_all_pairs(pairs):
        return [(score_offsets[i], score_offsets[j], fam) for i, j, fam in pairs]

    api_data = viewer_export.build_multichain_api_data(
        colored_json,
        result.chain_of,
        result.auth_of,
        structure_id,
        colored_svg_path=colored_svg if colored_svg.exists() else None,
    )
    # Model panel only: grey out nucleotides outside score_positions (e.g. a
    # widened reference's unmodeled second chain) using the existing
    # unobserved/unresolved-residue dimming path (backbone + text), so it
    # reads as "not part of this model" rather than looking like normal,
    # clickable model content. The reference panel keeps the plain api_data
    # — those residues are real reference data, nothing to greyed there.
    unscored_labels = sorted(p + 1 for p in range(n) if p not in score_positions)
    model_api_data = (
        {**api_data, "unobserved_label_seq_ids": unscored_labels}
        if unscored_labels
        else api_data
    )
    ref_fr3d = viewer_export.build_pairs_fr3d_data(
        result.nested_pairs,
        result.crossing_pairs,
        result.sequence,
        result.auth_of,
        structure_id,
        all_pairs=result.all_pairs,
    )

    # Full pairs (unscoped) for the reference panel's own display/list — shown
    # as-is even when widened, so the diagram honestly reflects the whole
    # structure. Scoped for the *scored* diff below, so a wider display never
    # changes the official matched/lost/added counts.
    ref_pairs = result.nested_pairs + result.crossing_pairs
    scoped_ref_pairs = [
        (i, j) for i, j in ref_pairs if i in score_positions and j in score_positions
    ]

    model_result = None
    model_cif = None
    if model_file:
        # Real predicted model: extract its base pairs and place them on the
        # reference coordinates (approach B).  Requires the same sequence in
        # the same chain order (co-indexed); we validate the length here.
        model_id = Path(model_file).stem
        model_cif = ensure_cif(Path(model_file), output_path / "downloads", quiet)
        chain_list = (
            [c.strip() for c in model_chains.split(",") if c.strip()]
            if model_chains
            else None
        )
        model_result = multichain.assemble(
            str(model_cif),
            str(output_path / "model_extraction"),
            chains=chain_list,
            auto_order=False,
            quiet=quiet,
        )
        if model_result is None:
            rprint("[red]Error: could not extract base pairs from the model[/red]")
            return
        if len(model_result.sequence) != score_n:
            rprint(
                f"[red]Model/reference length mismatch "
                f"({len(model_result.sequence)} vs {score_n}); the two must share "
                f"the same sequence in the same chain order. Aborting compare.[/red]"
            )
            return
        # Positioned onto result's label space (score_offsets is identity
        # unless the reference display was widened beyond score_chains).
        model_nested = _remap_pairs(model_result.nested_pairs)
        model_crossing = _remap_pairs(model_result.crossing_pairs)
        model_all_pairs = _remap_all_pairs(model_result.all_pairs)
        model_is_simulated = False
    else:
        model_id = f"{structure_id}_model"
        model_pairs = multichain.simulate_model_pairs(ref_pairs, n, seed=simulate_seed)
        model_nested, model_crossing = multichain.max_noncrossing(model_pairs, n)
        # The simulated model only perturbs canonical pairs.
        model_all_pairs = [(i, j, "cWW") for i, j in model_pairs]
        model_is_simulated = True

    # Interaction Network Fidelity of the model against the reference. Scoped
    # to score_chains (a no-op filter when unwidened) so a wider *display*
    # never inflates or deflates the score with pairs the model never
    # attempted (e.g. a second, unpredicted chain of a dimer).
    scoped_ref_all_pairs = [
        (i, j, fam)
        for i, j, fam in result.all_pairs
        if i in score_positions and j in score_positions
    ]
    inf_metrics = multichain.compute_inf(scoped_ref_all_pairs, model_all_pairs)

    # Shared-layout model panel is drawn on the reference 2D diagram, so
    # Nested filtering must use crossing relative to the *reference* nested
    # backbone — not the model's own 3D-topology nested/crossing split
    # (still kept in model_nested / model_crossing for the own-layout panel).
    ref_nested_set = {(min(i, j), max(i, j)) for i, j in result.nested_pairs}

    def _crosses_ref_nested(i: int, j: int) -> bool:
        a, b = (i, j) if i < j else (j, i)
        for c, d in ref_nested_set:
            if (a < c < b < d) or (c < a < d < b):
                return True
        return False

    shared_nested = []
    shared_crossing = []
    seen_pair = set()
    for i, j, _fam in model_all_pairs:
        a, b = (i, j) if i < j else (j, i)
        if (a, b) in seen_pair:
            continue
        seen_pair.add((a, b))
        if _crosses_ref_nested(a, b):
            shared_crossing.append((a, b))
        else:
            shared_nested.append((a, b))
    model_fr3d = viewer_export.build_pairs_fr3d_data(
        shared_nested,
        shared_crossing,
        result.sequence,
        result.auth_of,
        model_id,
        all_pairs=model_all_pairs,
    )

    viewer_dir = output_path / "viewer"
    viewer_dir.mkdir(exist_ok=True)
    structure_name = f"{structure_id}.{actual_format}"
    # Keep only the analysed RNA chain(s) in the 3D reference so a deposited
    # entry's extra chains (other RNA copies, ligands, whole proteins) don't
    # render as additional structures beside the one the model overlays.
    if not multichain.write_structure_chains(
        str(source_structure_path),
        list(result.order),
        str(viewer_dir / structure_name),
        quiet=quiet,
    ):
        shutil.copyfile(source_structure_path, viewer_dir / structure_name)
    copy_viewer_assets(viewer_dir)

    # Superimpose the real predicted model onto the reference so the 3D pane can
    # overlay both structures (approach C: pre-aligned coordinates — the pinned
    # pdbe-molstar bundle exposes no in-browser transform).  Reference and model
    # are co-indexed (approach B), giving an exact per-residue atom correspondence.
    overlays = []
    superpose_rmsd = None
    superpose_n = None
    if model_file and not model_is_simulated:
        # ref_index walks score_offsets (not 0..n-1) so a widened result still
        # only superposes on the chains the model actually corresponds to.
        ref_index = [
            (
                (result.chain_of[p], result.auth_of[p])
                if p < len(result.auth_of) and result.auth_of[p] is not None
                else None
            )
            for p in score_offsets
        ]
        model_index = [
            (
                (model_result.chain_of[i], model_result.auth_of[i])
                if i < len(model_result.auth_of) and model_result.auth_of[i] is not None
                else None
            )
            for i in range(score_n)
        ]
        aligned_name = f"{model_id}.aligned.cif"
        aligned_path = viewer_dir / aligned_name
        sp = multichain.superpose_model_onto_reference(
            str(source_structure_path),
            str(model_cif),
            ref_index,
            model_index,
            str(aligned_path),
            quiet=quiet,
        )
        if sp is None:
            # Too few correspondences to fit — show the model in its own frame.
            shutil.copyfile(str(model_cif), str(aligned_path))
        else:
            superpose_rmsd, superpose_n = sp
        # The model has its own author numbering, so the 3D selection needs a
        # label→(auth, chain) map built from the model (not the reference).
        # Keyed by both shared-layout labels (score_offsets[i]+1) and own-layout
        # labels (i+1) so either model 2D mode can drive the overlay.
        model_label_to_auth = {}
        model_label_to_chain = {}
        for i in range(score_n):
            shared_key = str(score_offsets[i] + 1)
            own_key = str(i + 1)
            if i < len(model_result.auth_of) and model_result.auth_of[i] is not None:
                model_label_to_auth[shared_key] = model_result.auth_of[i]
                model_label_to_auth[own_key] = model_result.auth_of[i]
            if i < len(model_result.chain_of):
                model_label_to_chain[shared_key] = model_result.chain_of[i]
                model_label_to_chain[own_key] = model_result.chain_of[i]
        # label-maps.json lives in model/ alongside that panel's own data (the
        # panels' ref/ and model/ dirs are created below, but this write runs
        # first, so make sure model/ exists here too).
        model_dir = viewer_dir / "model"
        model_dir.mkdir(exist_ok=True)
        (model_dir / "label-maps.json").write_text(
            json.dumps(
                {
                    "labelToAuth": model_label_to_auth,
                    "labelToChain": model_label_to_chain,
                }
            )
        )
        overlays.append(
            {
                "structureId": model_id,
                "structureUrl": aligned_name,
                "structureFormat": "cif",
                # Best predicted model: dark blue (CASP16 scheme).
                "baseColor": _CASP_MODEL_BLUE,
                "baseUrl": "model/",
            }
        )

    chains_label = "+".join(result.order)
    model_tag = "model, simulated" if model_is_simulated else "model"

    # The viewer uses structureId as the 2D plugin's pdbId, which is interpolated
    # into CSS selectors (e.g. querySelector('.rnaTopoSvg_<id>'),
    # '.rnaview_<id>_<n>') by both the plugin and the glue code. A dot in the id
    # (e.g. from a "<model>.trim.pdb" filename → stem "<model>.trim") makes the
    # selector parse the '.' as a *new class* (".rnaTopoSvg_x.trim" = class
    # rnaTopoSvg_x AND class trim), so those lookups match nothing and that panel's
    # selection/click path silently breaks. Sanitise the id to [A-Za-z0-9_].
    def _safe_id(value):
        return re.sub(r"[^A-Za-z0-9_]", "_", str(value))

    # Per-panel base-pair-list annotation: classify each listed pair against the
    # *other* structure's pair set (position keys in the shared label space).
    # Reference list → TP (also in model) / FN (missed); model list → TP / FP.
    def _pair_keys(pairs):
        return sorted({f"{min(i, j) + 1}_{max(i, j) + 1}" for i, j in pairs})

    # Compare over *all* families (cWW + non-canonical), by residue-pair
    # position, so the base-pair list's TP/FP/FN badges are correct for the
    # non-canonical pairs now shown too — not just the cWW backbone.
    ref_pair_keys = _pair_keys((i, j) for i, j, _ in result.all_pairs)
    model_pair_keys = _pair_keys((i, j) for i, j, _ in model_all_pairs)

    # Write each panel's data as sibling JSON files instead of inlining it into
    # index.html: the client already fetches api.json/fr3d.json/lbn.json from a
    # panel's own baseUrl when they aren't passed inline (resolvePanelData() in
    # r2dt-2d-3d-viewer.js), so this is a pure generation-side change. Keeps
    # index.html small and the raw data independently fetchable/cacheable.
    ref_dir = viewer_dir / "ref"
    model_dir = viewer_dir / "model"
    ref_dir.mkdir(exist_ok=True)
    model_dir.mkdir(exist_ok=True)
    ref_lbn = lbn_export.build_lbn_data(api_data, ref_fr3d)
    model_lbn = lbn_export.build_lbn_data(model_api_data, model_fr3d)
    (ref_dir / "api.json").write_text(json.dumps(api_data))
    (ref_dir / "fr3d.json").write_text(json.dumps(ref_fr3d))
    (ref_dir / "lbn.json").write_text(json.dumps(ref_lbn))
    (model_dir / "api.json").write_text(json.dumps(model_api_data))
    (model_dir / "fr3d.json").write_text(json.dumps(model_fr3d))
    (model_dir / "lbn.json").write_text(json.dumps(model_lbn))
    # Also drop a root-level copy of the reference's data (duplicated from
    # ref/): the shared Mol* pane below always links panelIndex 0 (reference)
    # and resolves its own api.json/fr3d.json relative to the viewer root, not
    # the panel's baseUrl, so those two files need to exist there too.
    (viewer_dir / "api.json").write_text(json.dumps(api_data))
    (viewer_dir / "fr3d.json").write_text(json.dumps(ref_fr3d))
    # otherPairKeys is the *other* panel's base-pair position keys, used only
    # to classify this panel's own list as TP/FP/FN -- same size concern as
    # apiData/fr3dData, so it goes in bp-compare.json alongside them rather
    # than inline in the panel dict.
    (ref_dir / "bp-compare.json").write_text(json.dumps(model_pair_keys))
    (model_dir / "bp-compare.json").write_text(json.dumps(ref_pair_keys))

    # Independent model layout (always for real models) so the viewer can
    # switch the right-hand 2D between shared-reference and own coordinates.
    model_own_layout = None
    if model_result is not None and not model_is_simulated:
        if not quiet:
            rprint("Generating model's own 2D layout (templatefree)...")
        own_layout = layout_multichain_structure(
            ctx,
            structure_id=_safe_id(model_id),
            sequence=model_result.sequence,
            dot_bracket=model_result.dot_bracket,
            boundaries=model_result.boundaries,
            nested_pairs=model_result.nested_pairs,
            crossing_pairs=model_result.crossing_pairs,
            layout_dir=output_path / "model_layout",
            rnapuzzler_flag=rnapuzzler_flag,
            quiet=quiet,
        )
        if own_layout is None:
            if not quiet:
                rprint(
                    "[yellow]Model own-layout generation failed; "
                    "shared layout only.[/yellow]"
                )
        else:
            own_json, own_svg = own_layout
            own_api = viewer_export.build_multichain_api_data(
                own_json,
                model_result.chain_of,
                model_result.auth_of,
                _safe_id(model_id),
                colored_svg_path=own_svg,
            )
            own_fr3d = viewer_export.build_pairs_fr3d_data(
                model_result.nested_pairs,
                model_result.crossing_pairs,
                model_result.sequence,
                model_result.auth_of,
                _safe_id(model_id),
                all_pairs=model_result.all_pairs,
            )
            # Ref pair keys remapped into the model's native 1..score_n space
            # so TP/FP/FN badges stay correct on the own-layout pair list.
            inv = {score_offsets[i]: i for i in range(score_n)}
            ref_keys_own = _pair_keys(
                (inv[i], inv[j])
                for i, j, _ in scoped_ref_all_pairs
                if i in inv and j in inv
            )
            own_dir = viewer_dir / "model-own"
            own_dir.mkdir(exist_ok=True)
            own_lbn = lbn_export.build_lbn_data(own_api, own_fr3d)
            (own_dir / "api.json").write_text(json.dumps(own_api))
            (own_dir / "fr3d.json").write_text(json.dumps(own_fr3d))
            (own_dir / "lbn.json").write_text(json.dumps(own_lbn))
            (own_dir / "bp-compare.json").write_text(json.dumps(ref_keys_own))
            # Dual-key maps (same as model/) so 3D click-through works.
            if (model_dir / "label-maps.json").is_file():
                shutil.copyfile(
                    model_dir / "label-maps.json", own_dir / "label-maps.json"
                )
            to_shared = {str(i + 1): score_offsets[i] + 1 for i in range(score_n)}
            from_shared = {str(score_offsets[i] + 1): i + 1 for i in range(score_n)}
            model_own_layout = {
                "baseUrl": "model-own/",
                "labelBridge": {"toShared": to_shared, "fromShared": from_shared},
            }
            if not quiet:
                rprint(f"[green]Model own layout ready: {own_svg}[/green]")

    panels = [
        {
            "title": f"{structure_id} (reference)",
            "subtitle": f"2D · chains {chains_label}",
            "structureId": _safe_id(structure_id),
            "chainId": "",
            "baseUrl": "ref/",
            "structureUrl": f"../{structure_name}",
            "structureFormat": actual_format,
            "bpCompare": {"role": "reference"},
        },
        {
            "title": f"{model_id} ({model_tag})",
            "subtitle": f"2D · chains {chains_label}",
            "structureId": _safe_id(model_id),
            "chainId": "",
            "baseUrl": "model/",
            "structureUrl": f"../{structure_name}",
            "structureFormat": actual_format,
            "bpCompare": {"role": "model"},
        },
    ]
    if model_own_layout:
        panels[1]["layoutModes"] = {
            "shared": {"baseUrl": "model/", "labelBridge": None},
            "own": model_own_layout,
        }
        panels[1]["defaultLayoutMode"] = "shared"
    if chain_views:
        panels[0]["chainViews"] = chain_views
    molstar = {
        "panelIndex": 0,
        "baseUrl": ".",
        "structureId": structure_id,
        "structureUrl": structure_name,
        "structureFormat": actual_format,
    }
    if overlays:
        # Experimental structure: RNA green (CASP16 scheme). These are RNA targets,
        # so a single RNA colour is used; DNA/ligand/protein colours are defined
        # above for when mixed-molecule targets get per-component colouring.
        molstar["baseColor"] = _CASP_RNA_GREEN
        molstar["overlays"] = overlays
        # Optional panel baseColor is structure-identity chrome only (3D /
        # legends). 2D letter selection stays orange so it does not collide
        # with TP green / FN blue on pair strokes.
        panels[0]["baseColor"] = _CASP_RNA_GREEN
        panels[1]["baseColor"] = _CASP_MODEL_BLUE
    heading = f"{structure_id} — reference vs {model_tag}"
    # Score-chain boundaries in the (possibly widened) reference label space.
    if score_chains:
        by_chain = {cid: (start, end) for cid, start, end in result.boundaries}
        score_boundaries = [
            (cid, by_chain[cid][0], by_chain[cid][1])
            for cid in score_chains
            if cid in by_chain
        ]
    else:
        score_boundaries = list(result.boundaries)
    inf_report = multichain.build_inf_report(
        structure_id=structure_id,
        model_id=model_id,
        chains=[cid for cid, _s, _e in score_boundaries],
        boundaries=score_boundaries,
        ref_pairs=scoped_ref_all_pairs,
        model_pairs=model_all_pairs,
        inf=inf_metrics,
        extra={
            "model_simulated": model_is_simulated,
            "model_own_layout": bool(model_own_layout),
        },
    )
    inf_scopes = inf_report.get("scopes") or []
    html_path = viewer_html.render_compare(
        viewer_dir,
        page_title=f"{structure_id} — reference vs model",
        heading=heading,
        subtitle=f"chains {chains_label} · shared 3D",
        panels=panels,
        molstar=molstar,
        metrics=inf_metrics,
        scopes=inf_scopes,
    )

    # Structured metrics for the batch dashboard (avoids parsing stdout).
    # Scoped, like INF above, so a widened display doesn't change the score.
    matched, lost, added = multichain.diff_pairs(
        scoped_ref_pairs, model_nested + model_crossing
    )

    def _inf_block(key):
        m = inf_metrics.get(key) or {}
        return {k: m.get(k) for k in ("inf", "ppv", "sty", "tp", "fp", "fn")}

    def _scope_summary(scope):
        inf = scope.get("inf") or {}
        keys = ("inf", "ppv", "sty", "tp", "fp", "fn")

        def _block(metric):
            m = inf.get(metric) or {}
            return {k: m.get(k) for k in keys}

        return {
            "id": scope.get("id"),
            "label": scope.get("label"),
            "type": scope.get("type"),
            "chains": scope.get("chains"),
            "inf": {k: _block(k) for k in ("wc", "nwc", "all")},
        }

    metrics_payload = {
        "structure_id": structure_id,
        "model_id": model_id,
        "model_simulated": model_is_simulated,
        "model_own_layout": bool(model_own_layout),
        "chains": [cid for cid, _s, _e in score_boundaries],
        "boundaries": [
            {"chain": cid, "start": start, "end": end}
            for cid, start, end in score_boundaries
        ],
        "inf": {k: _inf_block(k) for k in ("wc", "nwc", "all")},
        "scopes": [_scope_summary(s) for s in inf_scopes],
        "superpose_rmsd": superpose_rmsd,
        "superpose_n_atoms": superpose_n,
        "diff": {
            "matched": len(matched),
            "lost": len(lost),
            "added": len(added),
        },
    }
    (viewer_dir / "metrics.json").write_text(
        json.dumps(metrics_payload, indent=2) + "\n"
    )
    (viewer_dir / "inf-pairs.json").write_text(json.dumps(inf_report, indent=2) + "\n")
    (viewer_dir / "inf-pairs.csv").write_text(multichain.inf_report_to_csv(inf_report))

    if not quiet:
        rprint(
            f"[cyan]Base-pair diff (reference vs {model_id}): "
            f"{len(matched)} matched, {len(lost)} missing in model, "
            f"{len(added)} model-only[/cyan]"
        )

        def _fmt(metric):
            val = metric.get("inf")
            return "n/a" if val is None else f"{val:.3f}"

        rprint(
            "[cyan]INF (interaction network fidelity): "
            f"WC {_fmt(inf_metrics['wc'])}, "
            f"non-WC {_fmt(inf_metrics['nwc'])}, "
            f"all {_fmt(inf_metrics['all'])}[/cyan]"
        )
        rprint(f"[green]Compare viewer ready: {html_path.resolve()}[/green]")
        rprint(
            "[dim]Serve over HTTP, e.g.:\n"
            f"  python3 -m http.server -d {viewer_dir.resolve()} 8000\n"
            "then open http://localhost:8000/[/dim]"
        )


def copy_viewer_assets(viewer_dir: Path) -> None:
    """Copy the vendored viewer assets next to ``index.html``.

    The pdb-rna-viewer compiled bundle isn't on a CDN, isn't on npm, and
    the GitHub release downloads are served with ``application/octet-stream``
    which browsers refuse to load as a stylesheet. So we vendor it (plus
    the ``r2dt-2d-3d-viewer.js`` interaction glue) under ``data/viewer/`` in the R2DT
    repo (Apache-2.0) and copy it into each output folder.
    """
    src = Path(__file__).resolve().parents[1] / "data" / "viewer"
    wanted = (
        viewer_html.VIEWER_PLUGIN_FILENAME,
        viewer_html.VIEWER_CSS_FILENAME,
        viewer_html.R2DT_CSS_FILENAME,
        viewer_html.VIEWER_JS_FILENAME,
    )
    missing = [name for name in wanted if not (src / name).is_file()]
    if missing:
        raise click.ClickException(
            f"Missing vendored viewer assets in {src}: {', '.join(missing)}. "
            "The R2DT checkout looks incomplete."
        )
    # License texts travel with the assets: pdb-rna-viewer's Apache-2.0 (with
    # EMBL-EBI acknowledgement clause) and the tslib banner from its build.
    licenses = tuple(path.name for path in src.glob("*.LICENSE*"))
    for name in wanted + licenses:
        shutil.copyfile(src / name, viewer_dir / name)
