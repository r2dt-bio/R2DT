# Updating documentation

The R2DT documentation is generated using [Sphinx](https://www.sphinx-doc.org/en/master/) and [MyST](https://myst-parser.readthedocs.io/en/latest/).

To generate documentation locally:

1. Launch Docker or install requirements from `requirements.txt` (see [installation instructions](./installation.md)).

2. Generate html files:
    ```bash
    # from the root of the repository
    cd docs
    make html

    # or to automatically rebuild on changes:
    sphinx-autobuild docs docs/_build/html
    ```

3. Open `_build/index.html` in your browser.

To check that all URLs are live:
```bash
make linkcheck
```
