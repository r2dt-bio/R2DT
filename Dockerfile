ARG BASE_IMAGE_VERSION=v1.4

FROM rnacentral/r2dt-base:${BASE_IMAGE_VERSION}

# Create venv
ENV VENV=$RNA/venv
ENV PATH="$VENV/bin:$PATH"
RUN python3 -m venv $VENV && pip3 install --upgrade pip

ADD requirements.txt /tmp/requirements.txt
RUN pip3 install -r /tmp/requirements.txt --no-cache-dir

WORKDIR /rna/r2dt

ADD . /rna/r2dt
