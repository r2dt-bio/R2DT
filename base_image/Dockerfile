FROM gcc:13 as build

ENV tRNAscan_VER=2.0.11 \
    BioEasel_COMMIT=e7fae0ab43fc47058183b71ff498ee0d0d2de6a7 \
    BioEasel_VER=0.11 \
# Note: Changing version from 2.4.x to other versions may require changes to the ViennaRNA installation code below
    ViennaRNA_VER=2.6.3 \
    autoconf_VER=2.71 \
    RNA=/rna


FROM build AS rscape
# Install R-scape
RUN wget http://rivaslab.org/software/rscape/rscape.tar.gz && \
    tar -xzf rscape.tar.gz && \
    mv rscape_* rscape
RUN cd rscape && ./configure -q && make -j 4 && make install DESTDIR=/install/rscape


FROM build AS tRNAscan-SE
# Install tRNAScan-SE
RUN \
    wget https://github.com/UCSC-LoweLab/tRNAscan-SE/archive/v${tRNAscan_VER}.tar.gz && \
    tar -xzf v${tRNAscan_VER}.tar.gz && \
    rm v${tRNAscan_VER}.tar.gz
RUN cd tRNAscan-SE-${tRNAscan_VER} && ./configure -q && make && make install DESTDIR=/install/tRNAscan-SE


FROM build AS Bio-Easel
# Install Bio-Easel
RUN \
    cpan -T  Inline Inline::C
RUN \
    git clone https://github.com/nawrockie/Bio-Easel.git && \
    cd Bio-Easel && git checkout ${BioEasel_COMMIT} && \
    mkdir src && cd src && \
    wget -O easel-Bio-Easel.zip https://github.com/EddyRivasLab/easel/archive/Bio-Easel-${BioEasel_VER}.zip && \
    unzip easel-Bio-Easel.zip && \
    mv easel-Bio-Easel-${BioEasel_VER} easel
RUN \
    cd Bio-Easel && perl Makefile.PL
RUN \
    cd Bio-Easel && make &&  make install DESTDIR=/install/bio-easel

FROM build AS easel
# Install develop branch of easel
RUN git clone -b develop https://github.com/EddyRivasLab/easel
RUN cd easel && \
    autoconf && ./configure -q && make -j 4 && make install DESTDIR=/install/easel


FROM build AS traveler
# Install Traveler
RUN git clone https://github.com/cusbg/traveler
RUN cd traveler/src && git checkout 0ee67cbb9c5aaa2f98340065fd047c9a8feea53e && \
    make build && mkdir /install && cp ../bin/traveler /install/traveler


FROM build AS scripts
# Install jiffy infernal hmmer scripts
RUN \
    git clone https://github.com/nawrockie/jiffy-infernal-hmmer-scripts.git && \
    cd jiffy-infernal-hmmer-scripts && \
    git checkout 31d6c3b826d432a30620507830749cee58e15e68 && \
    chmod +x *.pl && \
    mkdir -p /install/jiffy-infernal-hmmer-scripts && \
    cp *.pl /install/jiffy-infernal-hmmer-scripts


FROM debian:bookworm-slim as final-build

ENV tRNAscan_VER=2.0.11 \
    BioEasel_COMMIT=e7fae0ab43fc47058183b71ff498ee0d0d2de6a7 \
    BioEasel_VER=0.09 \
# Note: Changing version from 2.4.x to other versions may require changes to the ViennaRNA installation code below
    ViennaRNA_VER=2.6.3 \
    autoconf_VER=2.71 \
    RNA=/rna

RUN apt-get -qq update && apt-get install -qq -y --no-install-recommends  \
    python3 python3-pip python3-venv \
    && rm -rf /var/lib/{apt,dpkg,cache,log}

RUN mkdir $RNA
WORKDIR $RNA

# Install Ribovore and Infernal
RUN \
    mkdir -p $RNA/ribovore && \
    cd $RNA/ribovore && \
    curl -sSL https://raw.githubusercontent.com/ncbi/ribovore/master/install.sh | bash -s -- "linux" && \
    cd $RNA/ribovore && \
    rm -Rf ncbi-blast vecscreen_plus_taxonomy rRNA_sensor

# Make sure tRNAScan-SE can find Infernal
RUN \
    ln -s /rna/ribovore/infernal/src/cmsearch /usr/local/bin/cmsearch && \
    ln -s /rna/ribovore/infernal/src/cmscan /usr/local/bin/cmscan

# Create venv
RUN python3 -m venv $RNA/venv && \
    . $RNA/venv/bin/activate && \
    pip3 install --upgrade pip && \
    pip3 install --upgrade viennarna

ENV \
    RIBOINSTALLDIR="$RNA/ribovore" \
    RIBOEASELDIR="$RNA/ribovore/infernal/easel/miniapps" \
    RIBOTIMEDIR=/usr/bin \
    RIBODIR="$RNA/ribovore/ribovore" \
    BIOEASELDIR="$RNA/Bio-Easel/blib/lib:$RNA/Bio-Easel/blib/arch:$RNA/Bio-Easel:$RNA/Bio-Easel/lib" \
    LC_ALL="C.UTF-8" LANG="C.UTF-8" \
    DATAPATH="/rna/RNAstructure/data_tables/"

ENV \
    RIBOSCRIPTSDIR="$RIBOINSTALLDIR/ribovore" \
    RIBOINFERNALDIR="$RIBOINSTALLDIR/infernal/src" \
    RIBOSEQUIPDIR="$RIBOINSTALLDIR/sequip" \
    RIBOBLASTDIR="$RIBOINSTALLDIR/ncbi-blast/bin" \
    VECPLUSDIR="$RIBOINSTALLDIR/vecscreen_plus_taxonomy" \
    BLASTDB="$VECPLUSDIR/univec-files":"$RRNASENSORDIR":"$BLASTDB" \
    RRNASENSORDIR="$RIBOINSTALLDIR/rRNA_sensor"

ENV \
    RIBOSCRIPTSDIR="$RIBOINSTALLDIR/ribovore" \
    RIBOINFERNALDIR="$RIBOINSTALLDIR/infernal/binaries" \
    RIBOEASELDIR="$RIBOINSTALLDIR/infernal/binaries" \
    RIBOSEQUIPDIR="$RIBOINSTALLDIR/sequip" \
    RIBOBLASTDIR="$RIBOINSTALLDIR/ncbi-blast/bin" \
    RIBOTIMEDIR=/usr/bin \
    RRNASENSORDIR="$RIBOINSTALLDIR/rRNA_sensor" \
    VECPLUSDIR="$RIBOINSTALLDIR/vecscreen_plus_taxonomy" \
    BLASTDB="$VECPLUSDIR/univec-files":"$RRNASENSORDIR":"$BLASTDB" \
    PERL5LIB="$RIBOSCRIPTSDIR":"$RIBOSEQUIPDIR":"$VECPLUSDIR":"$PERL5LIB" \
    PATH="$RIBOSCRIPTSDIR":"$RIBOBLASTDIR":"$RRNASENSORDIR":"$PATH"

ENV \
    PATH="$RIBOSCRIPTSDIR":"$RIBOBLASTDIR":"$RRNASENSORDIR":"$PATH" \
    PATH="$RNA/traveler/bin:$RIBODIR:$RIBOEASELDIR:$RIBOINFERNALDIR:/rna/rscape/bin:/rna/jiffy-infernal-hmmer-scripts:/rna/r2dt:/rna/RNAstructure/exe:$PATH" \
    PERL5LIB="$RIBOSCRIPTSDIR:$RIBOSEQUIPDIR:$VECPLUSDIR:$BIOEASELDIR:$RIBODIR:$EPNOPTDIR:$EPNOFILEDIR:$EPNTESTDIR:$PERL5LIB"

COPY --from=rscape /install/rscape /
COPY --from=tRNAscan-SE /install/tRNAscan-SE /
COPY --from=easel /install/easel /
COPY --from=Bio-Easel /install/bio-easel/* /
COPY --from=scripts /install/jiffy-infernal-hmmer-scripts/ ${RNA}/jiffy-infernal-hmmer-scripts/
COPY --from=traveler /install/traveler /usr/local/bin/traveler
