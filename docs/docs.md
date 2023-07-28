# Updating documentation

The R2DT documentation is generated using [Sphinx](https://www.sphinx-doc.org/en/master/) using [MyST](https://myst-parser.readthedocs.io/en/latest/).

To generate documentation locally:

1. Launch Docker or install requirements from `requirements.txt` (see [installation instructions](./installation.md)).

2. Generate html files:
    ```bash
    cd docs
    make html
    ```

3. Open `_build/index.html` in your browser.

To check that all URLs are live:
```bash
make linkcheck
```
