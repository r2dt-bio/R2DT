FROM rnacentral/r2dt-base:v1.4.1

# Create venv
ENV VENV=$RNA/venv
ENV PATH="$VENV/bin:$PATH"
RUN python3 -m venv $VENV && pip3 install --upgrade pip

ADD requirements.txt /tmp/requirements.txt
RUN pip3 install -r /tmp/requirements.txt --no-cache-dir

WORKDIR /rna/r2dt

ADD . /rna/r2dt
