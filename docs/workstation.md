# Local workstation

The R2DT **workstation** is a private local web UI for drawing sequences, structures, comparisons, and alignments. Jobs are cached on disk, can be edited (for structure modes), and exported either as a portable work package or as shareable static HTML.

It is complementary to the [command-line interface](./usage.md): the CLI remains the source of truth for pipeline options; the workstation is the interactive operator front-end.

## Requirements

- [Docker](https://www.docker.com) (or a compatible engine) with the `rnacentral/r2dt` image available
- A checkout of the R2DT repository (for `just workstation`)

Jobs run inside the same Docker image as the CLI, so Infernal, Traveler, FR3D, and the other R2DT tools do not need a separate host install.

## Starting the workstation

From the repository root:

```bash
just workstation
```

This mounts your workspace, publishes the UI on **http://127.0.0.1:8765/** (localhost only), and opens a browser when ready. Stop with `Ctrl+C`.

Equivalent CLI entry point (usually you do not need this if you use `just`):

```bash
r2dt.py workstation --workspace ~/.r2dt-workstation --port 8765
```

That defaults to `--bind 127.0.0.1`. Use `--bind 0.0.0.0` only inside Docker together with `-p 127.0.0.1:PORT:PORT` (as `just workstation` does) — never on a host that would listen on the LAN.

Environment overrides:

| Variable | Effect |
| --- | --- |
| `R2DT_WORKSPACE` | Cache directory (default `~/.r2dt-workstation`) |
| `R2DT_WORKSTATION_NO_OPEN=1` | Do not auto-open the browser |

## Trust model

The workstation is a **local operator tool**, not a multi-user web service:

- There is **no login or API token**. Anyone who can reach the HTTP port can create jobs (which run Docker), upload files, edit base pairs, and delete cached jobs.
- Prefer **`just workstation`**, which publishes the port to **127.0.0.1 only**. The process may bind `0.0.0.0` *inside* Docker so the published port works; that must not be confused with exposing the UI on the LAN.
- Mutating requests (`POST` / `PUT` / `DELETE`) require a loopback `Host` (`localhost` / `127.0.0.1` / `::1`). When a browser sends `Origin` or `Referer`, that must be loopback too. This blocks casual cross-origin CSRF and DNS rebinding; it is not a substitute for keeping the port off the network.
- Open the UI as `http://127.0.0.1:8765/` (or `http://localhost:8765/`). Stop the server when you are done.

Do not publish the port beyond loopback (for example `-p 0.0.0.0:8765:8765` or `--bind 0.0.0.0` on a host without Docker port filtering).

## Modes

The header switches between modes. Each mode has its own job list and “New …” form.

| Mode | URL | What it does |
| --- | --- | --- |
| **2D** | `/2d` | Draw secondary-structure diagrams from FASTA (one job per sequence for multi-FASTA) |
| **2D + 3D** | `/pdb` | Linked 2D topology + 3D Mol\* for one PDB/mmCIF |
| **2D + 2D + 3D** | `/compare` | Reference vs model comparison with INF scores, shared layout, and optional model own-layout |
| **Alignments** | `/align` | Stockholm / R-scape alignments → covariation-annotated diagrams |

Advanced options on each form mirror common CLI flags (layout mode, base-pair source, pseudoknots, and so on). Prefer the workstation form for interactive work; use the CLI docs when scripting.

## Jobs and cache

- Job metadata and outputs live under `~/.r2dt-workstation` (or `R2DT_WORKSPACE`).
- The dashboard lists jobs with **Status** and **Actions** next to the label; click the label to open a ready result.
- Import accepts a previously exported `.r2dt-job.zip`.
- CASP-style catalogs can be seeded into the cache with `just workstation-import-casp` when `output/site/casp15` / `casp16` already exist.

## Editing base pairs (structure modes)

On **2D + 3D** and **2D + 2D + 3D** viewers, the workstation can load an edit toolbar (add / delete / change Leontis–Westhof family). Edits are stored as overrides next to the job (`edits/`) and applied whenever you reopen that job. They are **not** written into the generated annotation files until you publish a shareable viewer (see below).

**2D** and **Alignments** results are static SVG galleries; they do not use the interactive pair editor.

## Nested toggle

In the 2D toolbar, **Nested** filters base pairs:

- **Off (default)** — show all pairs, including crossing / non-nested contacts
- **On** — show nested pairs only

On compare pages that use the **shared** (reference) layout, the model panel’s nested/crossing tags follow the reference backbone so the filter matches what you see on the diagram. The model’s own-layout view uses the model’s own nested/crossing split.

## Export

Ready jobs expose an **Export** menu with two targets:

1. **R2DT work package** (`.r2dt-job.zip`) — full job for another workstation: inputs, viewer or SVG results, edit overrides, and metadata. Re-import with **Import** on any mode dashboard.
2. **Shareable HTML** (`.r2dt-viewer.zip`) — static interactive viewer with edits **baked into** `fr3d.json` (and INF metrics refreshed when both panels exist). No workstation API required. Unzip, host on any static HTTPS origin, then open `index.html` or [embed](./pdb-2d-3d-viewer.md#embedding-on-your-own-site) with `R2DTViewer.create` / `createCompare`.

Shareable HTML is available only for jobs that have an interactive `viewer/` folder (typically **2D + 3D** and **2D + 2D + 3D**). Draw and alignment jobs export as work packages (and their SVG results pages).

## Related documentation

- [Command line reference](./usage.md) — `draw`, `pdb`, alignments, and other CLI commands
- [PDB structures](./pdb.md) — structure pipeline options
- [Interactive 2D + 3D viewer](./pdb-2d-3d-viewer.md) — `viewer/` layout, embed API, galleries
- [Stockholm alignments](./stockholm-alignments.md) — alignment inputs and R-scape
- [Docker images](./docker.md) — image build and upgrade notes
