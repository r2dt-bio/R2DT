# Updating documentation

The R2DT documentation is generated using [Sphinx](https://www.sphinx-doc.org/en/master/) and [MyST](https://myst-parser.readthedocs.io/en/latest/). The docs are deployed automatically to [Read the Docs](https://readthedocs.org/accounts/login/?next=/dashboard/).

To generate documentation locally:

1. Create a local virtual environment and activate it:
    ```bash
    just venv
    ```

2. Generate html files and automatically rebuild on changes:
    ```bash
    just docs
    ```

3. Open http://localhost:8000 in your browser to view the docs.

4. To check that all URLs are live:
    ```bash
    just check-links
    ```
