"""Multichain PDB pipeline helpers for the ``pdb`` / ``pdb_2d_3d``
commands: per-chain runs, combined layouts, RNAVIEW extraction.

Split out of r2dt.py mechanically; the CLI commands remain there.
"""

# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
# pylint: disable=too-many-branches,too-many-statements,duplicate-code

import json
import re
from pathlib import Path

from rich import print as rprint

from utils import compare_viewer, pdb_fetch
from utils import rnaview as rnaview_utils

# Set by r2dt.py after it defines the ``templatefree`` click command: these
# helpers re-run the 2D layout via ``ctx.invoke(templatefree, ...)``, and the
# command object lives in r2dt.py, which imports this module. Injecting the
# attribute avoids the circular import.
templatefree = None  # pylint: disable=invalid-name


def run_multichain_pdb(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    ctx,
    file_path,
    actual_format,
    structure_id,
    output_path,
    chains,
    rnapuzzler_flag,
    simulate_model,
    simulate_seed,
    compare,
    quiet,
    model_file=None,
    model_chains=None,
    _skip_partition=False,
):
    """Build a single combined 2D diagram from multiple RNA chains.

    Extracts the selected chains, concatenates them into one label space
    (auto-ordering unless an explicit order is given), lays out the nested pairs
    via the templatefree pipeline, then post-processes the SVG to break the
    phantom backbone bond at each chain junction, label per-chain 5'/3' termini,
    and draw crossing inter-chain pairs as an overlay.  The concatenation
    metadata is also written to ``<id>.multichain.json`` for the
    model-comparison draw.

    In ``--compare`` mode with an explicit ``--chains``, the requested chains
    may not be the whole story: the reference structure could have other RNA
    chains that interact with them (widen the display to the whole group —
    see ``multichain.partition_components``) or other chains entirely
    unrelated to them (offer a chain picker instead of silently dropping
    them). ``_skip_partition`` is for the sibling reference-only pages this
    generates for the picker — they already name an exact chain set and
    shouldn't recurse into building their own siblings.
    """
    # pylint: disable=import-outside-toplevel
    from utils import multichain

    if actual_format != "cif":
        rprint(
            "[red]Multi-chain mode requires an mmCIF structure "
            "(use --format cif or provide a .cif file).[/red]"
        )
        return

    if chains:
        chain_list = [c.strip() for c in chains.split(",") if c.strip()]
        auto_order = False  # explicit order provided
    else:  # --all-chains
        chain_list = None
        auto_order = True

    extraction_dir = output_path / "extraction"

    score_chains = None
    chain_views = None
    if compare and chain_list and not _skip_partition:
        components = multichain.partition_components(
            str(file_path), str(extraction_dir), quiet=True
        )
        if components:
            requested = set(chain_list)
            default_component = next(
                (c for c in components if requested & set(c)), None
            )
            if default_component:
                if set(default_component) != requested:
                    # The requested chains interact with others — widen the
                    # display to the whole group; scoring still only covers
                    # what was actually requested. This can happen even with
                    # only one component total (e.g. a 2-chain dimer where
                    # both chains interact and only one was requested).
                    score_chains = chain_list
                    chain_list = default_component
                    auto_order = True
                # A chain picker is only needed if there's something else in
                # the structure *not* part of the default view.
                other_components = [c for c in components if c != default_component]
                if other_components:
                    chain_views = build_chain_views(
                        ctx,
                        file_path,
                        actual_format,
                        structure_id,
                        output_path,
                        default_component,
                        other_components,
                        rnapuzzler_flag,
                        quiet,
                    )

    result = multichain.assemble(
        str(file_path),
        str(extraction_dir),
        chains=chain_list,
        auto_order=auto_order,
        quiet=quiet,
    )
    if result is None:
        rprint("[red]Error: multi-chain extraction failed[/red]")
        return

    counts = result.counts()
    if not quiet:
        rprint(
            f"Chains: {'+'.join(result.order)}  "
            f"({counts['length']} nt, {counts['nested']} nested pairs, "
            f"{counts['crossing']} crossing)"
        )

    # Record the concatenation metadata for later phases.
    meta_path = output_path / f"{structure_id}.multichain.json"
    with open(meta_path, "w") as meta_file:
        json.dump(
            {
                "structure_id": structure_id,
                "order": result.order,
                "boundaries": [
                    {"chain": c, "start": s, "end": e} for c, s, e in result.boundaries
                ],
                "sequence": result.sequence,
                "dot_bracket": result.dot_bracket,
                "nested_pairs": result.nested_pairs,
                "crossing_pairs": result.crossing_pairs,
                "counts": counts,
            },
            meta_file,
            indent=2,
        )

    # Write the combined FASTA (sequence + nested dot-bracket) and lay it out.
    fasta_path = output_path / f"{structure_id}.fasta"
    with open(fasta_path, "w") as fasta_file:
        fasta_file.write(f">{structure_id}\n{result.sequence}\n{result.dot_bracket}\n")

    results_folder = output_path / "results"
    if not quiet:
        rprint("Generating combined 2D diagram (templatefree)...")
    # Use auto mode (R2R + RNApuzzler + RNArtist, keep the fewest-overlap layout)
    # unless the caller forced RNApuzzler. Passing rscape=True would set chosen==1
    # and bypass the auto overlap-minimising selection, forcing R2R (which is not
    # overlap-free) -- so leave rscape False here.
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
    if svg_path is not None:
        # Break phantom junctions, add per-chain 5'/3' termini, draw the
        # crossing-pair overlay.
        multichain.postprocess_combined_svg(
            str(svg_path),
            result.boundaries,
            result.nested_pairs,
            result.crossing_pairs,
            quiet=quiet,
        )
        if not quiet:
            rprint(f"[green]Success! Combined SVG created: {svg_path}[/green]")

        if simulate_model and not compare:
            emit_simulated_model_panel(
                svg_path, structure_id, result, simulate_seed, quiet
            )
        if compare:
            compare_viewer.emit_compare_viewer(
                ctx,
                file_path,
                actual_format,
                structure_id,
                output_path,
                result,
                simulate_seed,
                quiet,
                model_file=model_file,
                model_chains=model_chains,
                score_chains=score_chains,
                chain_views=chain_views,
                rnapuzzler_flag=rnapuzzler_flag,
            )
    elif not quiet:
        rprint("[yellow]Diagram generation completed. Check output folder.[/yellow]")


def chain_group_slug(chain_ids):
    """Join chain ids into a filesystem-safe slug (e.g. ``A-B``)."""
    return re.sub(r"[^A-Za-z0-9]+", "-", "-".join(chain_ids)).strip("-") or "chain"


def build_chain_views(
    ctx,
    file_path,
    actual_format,
    structure_id,
    output_path,
    default_component,
    other_components,
    rnapuzzler_flag,
    quiet,
):
    """Generate a reference-only sibling page for each of ``other_components``
    — RNA chains of the same structure that don't interact with the ones
    shown by default — and return the ``chainViews`` list (the default page
    plus every sibling) for the reference panel's chain picker.

    Each sibling is a real (non-simulated) self-compare: the structure
    against itself, restricted to that component's chains. This reuses the
    existing compare-viewer machinery to get a real interactive 2D+3D page
    with no scoring semantics attached — INF=1.000/RMSD=0.00 there are the
    expected, meaningless artifact of comparing a structure to itself, not a
    real prediction score.
    """
    views = [
        {
            "label": f"Chain {'+'.join(default_component)} (current)",
            "url": "./index.html",
            "current": True,
        }
    ]
    for comp in other_components:
        slug = chain_group_slug(comp)
        sib_dir = output_path / f"chain-{slug}"
        sib_dir.mkdir(parents=True, exist_ok=True)
        run_multichain_pdb(
            ctx,
            file_path,
            actual_format,
            structure_id,
            sib_dir,
            ",".join(comp),
            rnapuzzler_flag,
            simulate_model=False,
            simulate_seed=2,
            compare=True,
            quiet=quiet,
            model_file=file_path,
            model_chains=",".join(comp),
            _skip_partition=True,
        )
        views.append(
            {
                # Sibling of the parent output dir, not of viewer/ itself — the
                # current page lives one level down, at <output>/viewer/index.html.
                "label": f"Chain {'+'.join(comp)}",
                "url": f"../chain-{slug}/viewer/index.html",
                "current": False,
            }
        )
    return views


def emit_simulated_model_panel(svg_path, structure_id, result, seed, quiet):
    """TESTING aid: draw a reference/model base-pair diff on the reference
    layout, using a randomly perturbed copy of the reference pairs as a stand-in
    model (see multichain.simulate_model_pairs)."""
    # pylint: disable=import-outside-toplevel
    from utils import multichain

    ref_pairs = result.nested_pairs + result.crossing_pairs
    model_pairs = multichain.simulate_model_pairs(
        ref_pairs, len(result.sequence), seed=seed
    )
    model_svg = svg_path.with_name(f"{structure_id}.model.svg")
    multichain.render_model_panel(
        str(svg_path), str(model_svg), ref_pairs, model_pairs, quiet=quiet
    )
    if not quiet:
        matched, lost, added = multichain.diff_pairs(ref_pairs, model_pairs)
        rprint(
            f"[cyan]Simulated model panel: {model_svg} "
            f"({len(matched)} matched, {len(lost)} missing, "
            f"{len(added)} model-only)[/cyan]"
        )


def rename_templated_outputs(results_folder: Path, structure_id: str):
    """Collapse ``<structure_id>-<template_id>.<ext>`` filenames produced
    by ``draw`` into plain ``<structure_id>.<ext>`` so the rest of the
    ``pdb`` / ``pdb_2d_3d`` pipeline (grey-out, viewer-export) doesn't
    need to know which template won.

    Returns the matched template id (the trailing portion after the
    structure id) if a templated SVG was found, otherwise ``None``.
    """
    svg_dir = results_folder / "results" / "svg"
    candidates = list(svg_dir.glob(f"{structure_id}-*.colored.svg"))
    if not candidates:
        return None
    # Take the first match -- there should be only one per structure.
    colored_svg = candidates[0]
    full_stem = colored_svg.name[: -len(".colored.svg")]
    template_id = full_stem[len(structure_id) + 1 :]

    suffix_dirs = {
        "results/svg": [".colored.svg", ".enriched.svg"],
        "results/thumbnail": [".thumbnail.svg"],
        "results/json": [".colored.json"],
        "results/fasta": [".fasta"],
    }
    for subdir, suffixes in suffix_dirs.items():
        d = results_folder / subdir
        if not d.is_dir():
            continue
        for suffix in suffixes:
            src = d / f"{full_stem}{suffix}"
            if src.exists():
                src.replace(d / f"{structure_id}{suffix}")
    return template_id


def extract_with_rnaview(pdb_file: str, chain_id=None, quiet=False):
    """
    Extract secondary structure using RNAView.

    Supports gzip-compressed PDB files (.pdb.gz).

    Args:
        pdb_file: Path to PDB file (may be gzip-compressed).
        chain_id: Optional chain ID. If None, uses first chain only.
        quiet: If True, suppress verbose output.

    Returns:
        Tuple of (sequence, dot_bracket) or (None, None).
    """
    try:
        # Use DecompressedStructureFile to handle .gz files
        # RNAView requires uncompressed file on disk
        with pdb_fetch.DecompressedStructureFile(Path(pdb_file)) as decompressed_path:
            # Extract sequence using rnaview module
            # If no chain specified, use first chain only (consistent with FR3D behavior)
            sequence = rnaview_utils.extract_sequence(
                str(decompressed_path), chain_id=chain_id, quiet=quiet
            )

            if not sequence:
                return None, None

            # Run RNAView
            rnaview_output = rnaview_utils.run_rnaview(str(decompressed_path))

            # Parse output to dot-bracket
            dot_bracket = rnaview_utils.parse_rnaview_output(
                rnaview_output, sequence, quiet
            )

            # Clean up temporary files created by RNAView
            rnaview_utils.cleanup_rnaview_files(str(decompressed_path))

            return sequence, dot_bracket

    except Exception as e:  # pylint: disable=broad-exception-caught
        if not quiet:
            print(f"Error in RNAView extraction: {e}")
        return None, None
