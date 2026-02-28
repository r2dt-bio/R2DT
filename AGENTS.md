# R2DT Copilot Instructions

R2DT visualizes RNA secondary structure using templates. It generates consistent, reproducible 2D diagrams by matching input sequences to covariance models from multiple RNA databases (CRW, Rfam, RiboVision, GtRNAdb, etc.).

## Architecture Overview

- **[r2dt.py](../r2dt.py)** - Main CLI entry point using Click. Defines all commands (`draw`, `pdb`, `rfam draw`, etc.). The `draw` command cascades through template libraries in priority order: RNAse P → tmRNA → RiboVision SSU/LSU → CRW → Rfam → GtRNAdb
- **[utils/core.py](../utils/core.py)** - Core `visualise()` function that runs the pipeline: cmalign → generate fasta/bpseq → call Traveler
- **[utils/config.py](../utils/config.py)** - All path constants for data directories and CM libraries
- **[utils/runner.py](../utils/runner.py)** - Subprocess executor wrapping shell commands. Use `runner.run(cmd)` for all external tool calls
- **[utils/shared.py](../utils/shared.py)** - Shared utilities including `get_ribotyper_output()` for template matching

## Development Workflow

**All development must use Docker** - the image includes Infernal, Traveler, RNAView, and other bioinformatics tools:

```bash
# Run tests (preferred method)
just test-all                    # Run full test suite
just test Draw                   # Run specific test class (e.g., TestDraw)

# Interactive shell for debugging
just run                         # docker run with mounted cwd
docker run --rm -v $(pwd):/rna/r2dt -it rnacentral/r2dt:latest bash

# Run commands inside container
python3 r2dt.py draw examples/examples.fasta output/ --quiet
python3 r2dt.py test             # Alternative to just test-all
```

## Testing Conventions

Tests in [tests/tests.py](../tests/tests.py) follow a pattern:
- Extend `R2dtTestCase` base class
- Define `cmd`, `test_results`, `precomputed_results`, `files` class attributes
- Override `test_examples()` to call `self.check_examples()` for SVG comparison
- Tests compare generated SVGs against reference images using SSIM (structural similarity)
- Set `R2DT_KEEP_TEST_RESULTS=1` to preserve test output for inspection

## Code Style

**All code must pass Black and Pylint before committing.**

- **Formatting**: Black (line-length 88), isort with "black" profile
- **Linting**: Pylint - use `# pylint: disable=...` comments for legitimate suppressions
- **Imports**: Standard library → Third-party (click, rich, requests) → Local (`from utils import ...`)
- **CLI output**: Use `rich.print` (`rprint`) for colored terminal output
- **Quiet mode**: Most commands have `--quiet/-q` flag; check before printing

```bash
# Run inside Docker container
black .                          # Format all files
black --check .                  # Check without modifying
pylint r2dt.py utils/ tests/     # Run linter
```

All code must score 10.0 on Pylint and be formatted by Black.

## Key Patterns

### Adding new RNA template sources
1. Add constants to [utils/config.py](../utils/config.py) for CM library paths
2. Create module in `utils/` (see [utils/rfam.py](../utils/rfam.py), [utils/gtrnadb.py](../utils/gtrnadb.py))
3. Add command group in [r2dt.py](../r2dt.py) using `@cli.group()` / `@group.command()`
4. Integrate into cascade in `draw()` command

### Working with covariance models
```python
from utils.shared import cmfetch, get_ribotyper_output
cm_file = cmfetch(model_id, cm_library)  # Fetches single CM from combined all.cm
```

### External tool execution
```python
from utils.runner import runner
runner.run(f"cmalign --outformat stockholm {cm_file} {fasta}")
```

## Data Directory Structure

```
data/
├── crw/           # CRW rRNA templates (cms)
├── rfam/          # Rfam families with traveler templates
│   └── {RF*}/     # Per-family: traveler-template.xml, *-traveler.fasta
├── gtrnadb/       # tRNA templates by domain/isotype
├── ribovision-lsu/# Large subunit templates
├── ribovision-ssu/# Small subunit templates
├── rnasep/        # RNase P templates
└── tmrna/         # tmRNA templates
```
