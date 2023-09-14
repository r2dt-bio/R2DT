FROM gcc:13 as build

ENV tRNAscan_VER=2.0.11 \
    BioEasel_COMMIT=e7fae0ab43fc47058183b71ff498ee0d0d2de6a7 \
    traveler_COMMIT=0ee67cbb9c5aaa2f98340065fd047c9a8feea53e \
    jiffy_COMMIT=31d6c3b826d432a30620507830749cee58e15e68 \
    BioEasel_VER=0.11


FROM build AS rscape

# Install R-scape
RUN \
    wget --no-verbose -O - http://rivaslab.org/software/rscape/rscape.tar.gz | tar zxf - 2>/dev/null && \
    mv rscape_* rscape
RUN cd rscape && ./configure -q && make -j 2 && make install DESTDIR=/install/


FROM build AS tRNAscan-SE

# Install tRNAScan-SE
RUN \
    wget --no-verbose -O - https://github.com/UCSC-LoweLab/tRNAscan-SE/archive/v${tRNAscan_VER}.tar.gz | tar zxf -
RUN \
    cd tRNAscan-SE-${tRNAscan_VER} && ./configure -q --prefix=/install/ && make && make install


FROM build AS Bio-Easel

# Install Bio-Easel
RUN \
    cpan -T  Inline Inline::C > /dev/null
RUN \
    git clone https://github.com/nawrockie/Bio-Easel.git && \
    cd Bio-Easel && git checkout ${BioEasel_COMMIT} && \
    mkdir src && cd src && \
    wget --no-verbose -O - https://github.com/EddyRivasLab/easel/archive/Bio-Easel-${BioEasel_VER}.tar.gz | tar zxf - && \
    mv easel-Bio-Easel-${BioEasel_VER} easel
RUN \
    cd Bio-Easel && perl Makefile.PL
RUN \
    cd Bio-Easel && make && make install DESTDIR=/install/ && \
    cp -rv /Bio-Easel/lib /install/Bio-Easel/lib && \
    cp -rv /Bio-Easel/blib/lib /install/Bio-Easel/blib/lib

FROM build AS traveler

# Install Traveler
RUN \
    git clone https://github.com/cusbg/traveler
RUN \
    cd traveler/src && git checkout ${traveler_COMMIT} && \
    make build && mkdir -p /install/utils && cp ../bin/traveler /install/traveler
RUN \
    cp -rv traveler/utils/* /install/utils/


FROM build AS scripts

# Install jiffy infernal hmmer scripts
RUN \
    git clone https://github.com/nawrockie/jiffy-infernal-hmmer-scripts.git && \
    cd jiffy-infernal-hmmer-scripts && git checkout ${jiffy_COMMIT} && \
    chmod +x *.pl && mkdir -p /install/ && cp *.pl /install/


FROM build as ribovore-infernal-easel

# Install Ribovore, Infernal, and easel
RUN \
    mkdir -p /install/ && \
    cd /install/ && \
    wget --no-verbose -O - https://raw.githubusercontent.com/ncbi/ribovore/install-r2dt/install.sh | bash -s -- "rt-linux" && \
    rm -Rf /install/ribovore/models/*.stk /install/ribovore/models/*.cm* /install/ribovore/testfiles

FROM debian:bookworm-slim as final-build

ENV \
    RNA=/rna

RUN \
    apt-get -qq update && apt-get install -qq -y --no-install-recommends  \
    python3  \
    python3-pip  \
    python3-venv  \
    python3-cairo  \
    perl &&  \
    rm -rf /var/lib/{apt,lpkg,cache,log}
RUN \
    cpan -T Time::HiRes

RUN mkdir $RNA
WORKDIR $RNA

## Make sure tRNAScan-SE can find Infernal
RUN \
    ln -s /rna/ribovore/infernal/binaries/cmsearch /usr/bin/cmsearch && \
    ln -s /rna/ribovore/infernal/binaries/cmscan /usr/bin/cmscan

ENV \
    RIBOINSTALLDIR="$RNA/ribovore" \
    RIBOEASELDIR="$RNA/ribovore/infernal/binaries" \
    RIBOTIMEDIR=/usr/bin \
    RIBODIR="$RNA/ribovore/ribovore" \
    BIOEASELDIR="$RNA/Bio-Easel/blib/lib:$RNA/Bio-Easel/blib/arch:$RNA/Bio-Easel:$RNA/Bio-Easel/lib" \
    LC_ALL="C.UTF-8" LANG="C.UTF-8" \
    DATAPATH="/rna/RNAstructure/data_tables/"

ENV \
    RIBOSCRIPTSDIR="$RIBOINSTALLDIR/ribovore" \
    RIBOINFERNALDIR="$RIBOINSTALLDIR/infernal/binaries" \
    RIBOSEQUIPDIR="$RIBOINSTALLDIR/sequip" \
    RIBOBLASTDIR="$RIBOINSTALLDIR/ncbi-blast/bin" \
    VECPLUSDIR="$RIBOINSTALLDIR/vecscreen_plus_taxonomy" \
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
    PERL5LIB="$RIBOSCRIPTSDIR":"$RIBOSEQUIPDIR":"$VECPLUSDIR":"$PERL5LIB" \
    PATH="$RIBOSCRIPTSDIR":"$RIBOBLASTDIR":"$RRNASENSORDIR":"$PATH"

ENV \
    PATH="$RIBOSCRIPTSDIR":"$RIBOBLASTDIR":"$RRNASENSORDIR":"$PATH" \
    PATH="$RNA/traveler/bin:$RIBODIR:$RIBOEASELDIR:$RIBOINFERNALDIR:/rna/rscape/bin:/rna/jiffy-infernal-hmmer-scripts:/rna/r2dt:/rna/RNAstructure/exe:$PATH" \
    PERL5LIB="$RIBOSCRIPTSDIR:$RIBOSEQUIPDIR:$VECPLUSDIR:$BIOEASELDIR:$RIBODIR:$EPNOPTDIR:$EPNOFILEDIR:$EPNTESTDIR:$PERL5LIB:/usr/lib/tRNAscan-SE"

COPY --from=rscape /install/ /{RNA}/rscape

COPY --from=tRNAscan-SE /install/ /usr
RUN sed -i -e 's|/install/bin|/usr/bin|' -e 's|/install/lib|/usr/lib|' /usr/bin/tRNAscan-SE.conf

COPY --from=Bio-Easel /install/usr/* /usr
COPY --from=Bio-Easel /install/Bio-Easel $RNA/Bio-Easel
COPY --from=Bio-Easel /usr/local/share/perl/5.36.0/ /usr/share/perl/5.36/

COPY --from=scripts /install/ ${RNA}/jiffy-infernal-hmmer-scripts/

COPY --from=ribovore-infernal-easel /install/ ${RNA}/ribovore

COPY --from=traveler /install/traveler /usr/local/bin/traveler
COPY --from=traveler /install/utils/* ${RNA}/traveler/utils/
# Forgive me, for I have hacked
# In the name of the Compiler, the Debugger, and the Holy Stack Trace
COPY --from=traveler "/usr/local/lib64/libstdc++.so.6.0.32" "/usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.30"
