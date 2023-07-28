# API

The [R2DT API](https://www.ebi.ac.uk/Tools/common/tools/help/index.html?tool=r2dt) is powered by the EMBL-EBI Web Services. The API allows submitting a job, checking the job status until it is complete, and then fetching the results in various formats.

For example, it is possible to download results in the [RNA 2D JSON Schema](https://github.com/LDWLab/RNA2D-data-schema/) format:

1. Submit a sequence to the [R2DT API](https://www.ebi.ac.uk/Tools/common/tools/help/index.html?tool=r2dt) and get a job id like `r2dt-R20230726-220347-0839-95483427-p1m`.

1. Once the results are ready, the JSON schema file can be fetched from a URL like this: `https://www.ebi.ac.uk/Tools/services/rest/r2dt/result/r2dt-R20230726-220347-0839-95483427-p1m/json`.

Note that the job ids expire within 7 days after the job is completed.
