# Updating documentation

The R2DT documentation is generated using [Sphinx](https://www.sphinx-doc.org/en/master/) and [MyST](https://myst-parser.readthedocs.io/en/latest/). The docs are deployed automatically to [Read the Docs](https://readthedocs.org/accounts/login/?next=/dashboard/).

To generate documentation locally:

1. Create a local virtual environment and activate it:
    ```bash
    python3 -m venv .venv
    pip3 install -r requirements.txt
    source .venv/bin/activate
    ```

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
cd docs
make linkcheck
```
