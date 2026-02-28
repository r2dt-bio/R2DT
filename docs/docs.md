# Updating documentation

The R2DT documentation is generated using [Sphinx](https://www.sphinx-doc.org/en/master/) and [MyST](https://myst-parser.readthedocs.io/en/latest/). The docs are deployed automatically to [Read the Docs](https://readthedocs.org/accounts/login/?next=/dashboard/).

To generate documentation locally:

1. Start local server and watch for changes:
    ```bash
    just docs
    ```

2. Open [http://0.0.0.0:8000](http://0.0.0.0:8000) in your browser to view the docs.

:::{note}
An alternative port can be specified using the `port` variable: `just port=3000 docs`
:::

## Regenerating doc images

The viral genome diagrams in `docs/images/` are generated from example inputs in `examples/viral/` and `examples/hcv-alignment.stk`. To regenerate all of them:

```bash
just docs-images
```

This runs `viral-annotate`, `stitch`, and `stockholm` inside Docker against the checked-in example FASTA and Stockholm files, then copies the resulting SVGs into `docs/images/`. Run this after any algorithm change that affects diagram output, and commit the updated images alongside the code change.

:::{tip}
Use `just docs-images some-tag` to regenerate against a specific Docker image tag.
:::

## Verifying documentation

To check that all hyperlinks are live:

    just check-links
