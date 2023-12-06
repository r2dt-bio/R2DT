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

## Verifying documentation

To check that all hyperlinks are live:

    just check-links
