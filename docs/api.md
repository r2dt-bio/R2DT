# API and Python package

## API

The [R2DT API](https://www.ebi.ac.uk/Tools/common/tools/help/index.html?tool=r2dt) is powered by the EMBL-EBI Web Services. The API allows submitting a job, checking the job status until it is complete, and then fetching the results in various formats.

For example, it is possible to download results in the [RNA 2D JSON Schema](https://github.com/LDWLab/RNA2D-data-schema/) format:

1. Submit a sequence to the [R2DT API](https://www.ebi.ac.uk/Tools/common/tools/help/index.html?tool=r2dt) and get a job id like `r2dt-R20230726-220347-0839-95483427-p1m`.

1. Once the results are ready, the JSON schema file can be fetched from a URL like this: `https://www.ebi.ac.uk/Tools/services/rest/r2dt/result/r2dt-R20230726-220347-0839-95483427-p1m/json`.

Note that the job ids expire within 7 days after the job is completed.

## Python package

The [r2dt-client](https://github.com/anayden/r2dt-client) Python package provides a convenient way to access the R2DT API from Python scripts, Jupiter notebooks, and Google Colab.

### Running r2dt-client in Google Colab

```
!pip install r2dt_client

from r2dt_client import setup, draw

setup(email="antonipetrov@gmail.com")
draw(
    ">S box leader\nCTCTTATCGAGAGTTGGGCGAGGGATTTGGCCTTTTGACCCCAAAAGCAACCGACCGTAATTCCATTGTGAAATGGGGCGCATTTTTTTCGCGCCGAGACGCTGGTCTCTTAAGGCACGGTGCTAATTCCATTCAGATCTGATCTGAGAGATAAGAG")
```
