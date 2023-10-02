# Testing

R2DT comes with a comprehensive testing suite (see [tests.py](https://github.com/RNAcentral/R2DT/blob/main/tests/tests.py)) which runs a predefined [set of sequences](https://github.com/RNAcentral/R2DT/tree/main/examples) through R2DT and compares the resulting SVG files with precomputed examples.

**All tests must pass** before a pull request can be merged into the `main` branch.

## Useful commands

* Run all tests
    ```bash
    r2dt.py test
    ```

    or

    ```bash
    python3 -m unittest
    ```

* Run a [single test](https://github.com/RNAcentral/R2DT/blob/main/tests/tests.py)
    ```bash
    r2dt.py test TestRibovisionLSU
    python3 -m unittest tests.tests.TestRibovisionLSU
    ```

* When running tests with `r2dt.py test` the results and intermediate files are kept by default to facilitate debugging. If tests are launched with `python3 -m unittest`, the results are deleted at the end of the tests. To keep the results, set the following environment variable:

    ```bash
    R2DT_KEEP_TEST_RESULTS=1 python3 -m unittest tests.tests.TestRnasep
    ```

## Updating example files

Sometimes a change in the algorithm results in a different secondary structure layout, and the tests are designed to alert the developers about any changes. However, sometimes changes are normal and expected, which means that precomputed examples need to be reviewed and updated, if needed.

When a test produces a structure that is not identical to the precomputed example, an **HTML report** including both structures is automatically generated and placed in the `tests` folder. The HTML files are named according to the name of the test.

For example, a file `TestRnasep_2.html` is generated when the second sequence in the `TestRnasep` class in [tests.py](https://github.com/RNAcentral/R2DT/blob/main/tests/tests.py) is different. These files are meant to be examined visually, and if the changes are acceptable, the examples can be updated with a helper command:

```bash
r2dt.py update-test-examples <test-name>
r2dt.py update-test-examples TestRfam
```

Once the precomputed examples are updated, it is recommended to check the output of `git diff` for the SVGs before committing the new example files, in case some changes were not picked up by eye.
