ARG BASE_IMAGE_VERSION=v.2.2.1

FROM rnacentral/r2dt-base:${BASE_IMAGE_VERSION}

WORKDIR /rna/r2dt

# Set up the Rfam and CRW data directories
ADD data/rfam/cms/all.cm.tar.gz data/rfam/cms/
RUN \
    rm -f data/rfam/cms/all.cm.ssi && \
    cmfetch --index data/rfam/cms/all.cm

ADD data/crw/all.cm.tar.gz data/crw/
RUN \
    rm -f data/crw/all.cm.ssi && \
    cmfetch --index data/crw/all.cm

# Create venv
ENV VENV=$RNA/venv
ENV PATH="$VENV/bin:$PATH"
RUN python3 -m venv $VENV && pip3 install --upgrade pip

ADD requirements.txt /tmp/requirements.txt
RUN pip3 install -r /tmp/requirements.txt --no-cache-dir

# Install FR3D-python
RUN pip3 install --no-cache-dir https://github.com/BGSU-RNA/fr3d-python/archive/3c7dc20.tar.gz

ADD . /rna/r2dt

# Index covariance model libraries included as plain files
RUN rm -f data/rnasep/cms/all.cm.ssi data/tmrna/cm/all.cm.ssi \
        data/ribovision-ssu/cms/all.cm.ssi data/ribovision-lsu/cms/all.cm.ssi && \
    cmfetch --index data/rnasep/cms/all.cm && \
    cmfetch --index data/tmrna/cm/all.cm && \
    cmfetch --index data/ribovision-ssu/cms/all.cm && \
    cmfetch --index data/ribovision-lsu/cms/all.cm
