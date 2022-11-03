# Updating documentation

The R2DT documentation is generated using [Sphinx](https://www.sphinx-doc.org/en/master/) using [MyST](https://myst-parser.readthedocs.io/en/latest/).

To generate documentation locally:
```bash
cd docs
make html
```
and open file `_build/index.html` in your browser.

To check that all URLs are live:
```bash
make linkcheck
```
