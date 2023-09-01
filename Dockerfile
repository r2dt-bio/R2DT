FROM rnacentral/r2dt-base

ADD . /rna/r2dt
RUN . $RNA/venv/bin/activate && pip3 install -r $RNA/r2dt/requirements.txt
