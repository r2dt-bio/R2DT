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

# Install FR3D-python, pinned to the tip of the actively-developed `latest`
# branch (bump the commit deliberately). That branch ships an fr3d/modified/
# subpackage but omits its __init__.py, so setup.py's find_packages() drops it
# and `import fr3d.cif.reader` fails (it pulls in fr3d.modified.mapping). Add the
# missing __init__.py before installing. The final import is a build-time guard.
ARG FR3D_COMMIT=ed850c00df01616e58c643b0f84bdc69662b5d55
RUN set -eux; \
    cd /tmp; \
    python3 -c "import urllib.request; urllib.request.urlretrieve('https://github.com/BGSU-RNA/fr3d-python/archive/${FR3D_COMMIT}.tar.gz', 'fr3d.tar.gz')"; \
    tar xzf fr3d.tar.gz; \
    touch "fr3d-python-${FR3D_COMMIT}/fr3d/modified/__init__.py"; \
    pip3 install --no-cache-dir "./fr3d-python-${FR3D_COMMIT}"; \
    rm -rf /tmp/fr3d.tar.gz "/tmp/fr3d-python-${FR3D_COMMIT}"; \
    python3 -c "from fr3d.cif.reader import Cif"

ADD . /rna/r2dt

# Index covariance model libraries included as plain files
RUN rm -f data/rnasep/cms/all.cm.ssi data/tmrna/cm/all.cm.ssi \
        data/ribovision-ssu/cms/all.cm.ssi data/ribovision-lsu/cms/all.cm.ssi && \
    cmfetch --index data/rnasep/cms/all.cm && \
    cmfetch --index data/tmrna/cm/all.cm && \
    cmfetch --index data/ribovision-ssu/cms/all.cm && \
    cmfetch --index data/ribovision-lsu/cms/all.cm
