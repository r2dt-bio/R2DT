# Vendored layered base-pair notation (LBN)

R2DT vendors a **minimal slice** of the standalone [layered base-pair
notation](https://github.com/na-hackathon/na-hackathon-2026/tree/main/workstreams/ws2-prediction-non-Watson-Crick/layered-bp-notation)
project here. That upstream repo is the full tool (CIF I/O, round-trip
examples, CLI); it is not on PyPI yet, so we copy the files R2DT needs
for [`--basepairs cif`](../cif_basepairs.py) rather than adding a dependency.

## Files in this directory

| File | Role |
| --- | --- |
| `common.py` | LW-family grouping, bracket packing, CIF-oriented helpers |
| `standalone_lbn_script.py` | mmCIF reader used by `utils/cif_basepairs.py` |

Only these two files are synced from upstream. Do not edit them in place for
R2DT-specific behaviour — change upstream or wrap from `utils/cif_basepairs.py`.

## License

The upstream repository is MIT-licensed (code) by the NA Hackathon 2026
contributors; the upstream `LICENSE` is vendored verbatim in this directory
alongside the code it covers. MIT code may be redistributed inside
Apache-2.0-licensed R2DT provided the notice is retained — which this
directory does. (The upstream `LICENSE-docs` CC BY 4.0 applies only to their
documentation, none of which is vendored here.)

## Refreshing from upstream

```bash
just sync-lbn
```

This downloads the latest `common.py` and `standalone_lbn_script.py` from the
hackathon repo path above.

## Related R2DT code (not vendored here)

- **`utils/cif_basepairs.py`** — wires the vendored script into the `pdb` /
  `pdb_2d_3d` pipeline when `--basepairs cif` is used.
- **`utils/lbn_export.py`** — builds `lbn.json` for the interactive viewer
  panel. Works in 1-based label space from `api.json` / `fr3d.json`; it
  reuses the same bracket-packing *idea* but is maintained separately because
  the viewer does not re-read mmCIF coordinates.

For the full LBN specification, examples, and `notation_to_cif.py`, see the
upstream project README in the hackathon repository.
