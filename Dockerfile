FROM rnacentral/r2dt-base

ADD . /rna/r2dt
ADD requirements.txt $RNA/r2dt/requirements.txt
RUN pip3 install -r $RNA/r2dt/requirements.txt
