FROM gcc:13 AS build

ENV tRNAscan_VER=2.0.11 \
    BioEasel_COMMIT=a52f0f7f7cae879bb95489bb794c04c1a0110aa7 \
    traveler_COMMIT=0ee67cbb9c5aaa2f98340065fd047c9a8feea53e \
    jiffy_COMMIT=31d6c3b826d432a30620507830749cee58e15e68 \
    Rscape_VER=v2.0.4.c \
    Ribovore_COMMIT=5df943c6f2866d7979bb13c948ca9d8531076898 \
    BioEasel_VER=0.11 \
    INFERNAL_VER=1.1.5 \
    RNArtist_VER=0.4.6


FROM build AS rscape

# Install R-scape
RUN \
    wget --no-verbose -O - http://rivaslab.org/software/rscape/rscape_${Rscape_VER}.tar.gz | tar zxf - 2>/dev/null && \
    mv rscape_${Rscape_VER} rscape
RUN  cp /usr/share/automake-1.16/config.guess /rscape/lib/R2R/R2R-current/config.guess && \
      mkdir -p /rna/rscape && cd rscape && ./configure --prefix=/rna/rscape -q && make -j 2 && make install


FROM build AS trnascan-se

# Install tRNAScan-SE
RUN \
    wget --no-verbose -O - https://github.com/UCSC-LoweLab/tRNAscan-SE/archive/v${tRNAscan_VER}.tar.gz | tar zxf -
RUN \
    cd tRNAscan-SE-${tRNAscan_VER} && ./configure -q --prefix=/install/ && make && make install


FROM build AS bio-easel

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


FROM build AS ribovore-infernal-easel

# Install Ribovore
RUN \
    mkdir -p /install/ && \
    cd /install/ && \
    git clone https://github.com/ncbi/ribovore && \
    cd ribovore && \
    git checkout ${Ribovore_COMMIT} && \
    rm -Rf .git models testfiles taxonomy

# Install Sequip (Ribovore dependency)
RUN \
    cd /install/ && \
    sequip_version=ribovore-1.0.5 && \
    wget -O sequip-${sequip_version}.zip https://github.com/nawrockie/sequip/archive/${sequip_version}.zip && \
    unzip sequip-${sequip_version}.zip && \
    mv sequip-${sequip_version} sequip && \
    rm sequip-${sequip_version}.zip

# Install Infernal
RUN \
    wget -q "http://eddylab.org/infernal/infernal-${INFERNAL_VER}.tar.gz" && \
    tar xzf "infernal-${INFERNAL_VER}.tar.gz" && \
    cd "infernal-${INFERNAL_VER}" && \
    mkdir /infernal && \
    ./configure --prefix=/infernal && \
    make && \
    make install && \
    cd easel && \
    make install && \
    mkdir -p /install/infernal/binaries && \
    mv /infernal/bin/* /install/infernal/binaries/

FROM build AS rnartist

# Install RNArtist
RUN \
    wget -q -O rnartist.jar https://github.com/fjossinet/RNArtistCore/releases/download/${RNArtist_VER}-SNAPSHOT/rnartistcore-${RNArtist_VER}-SNAPSHOT-jar-with-dependencies.jar

FROM debian:bookworm-slim AS final-build

ENV \
    RNA=/rna

RUN \
    apt-get -qq update && apt-get install -qq -y --no-install-recommends  \
    ncbi-blast+ \
    openjdk-17-jre-headless \
    python3  \
    python3-pip  \
    python3-venv  \
    python3-cairo  \
    chafa \
    procps \
    perl &&  \
    rm -rf /var/lib/{apt,lpkg,cache,log}
RUN \
    cpan -T Time::HiRes

RUN mkdir $RNA
WORKDIR $RNA

# Make sure tRNAScan-SE can find Infernal
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
    RIBOBLASTDIR="/usr/bin" \
    VECPLUSDIR="$RIBOINSTALLDIR/vecscreen_plus_taxonomy" \
    RRNASENSORDIR="$RIBOINSTALLDIR/rRNA_sensor"

ENV \
    RIBOSCRIPTSDIR="$RIBOINSTALLDIR/ribovore" \
    RIBOINFERNALDIR="$RIBOINSTALLDIR/infernal/binaries" \
    RIBOEASELDIR="$RIBOINSTALLDIR/infernal/binaries" \
    RIBOSEQUIPDIR="$RIBOINSTALLDIR/sequip" \
    RIBOTIMEDIR=/usr/bin \
    RRNASENSORDIR="$RIBOINSTALLDIR/rRNA_sensor" \
    VECPLUSDIR="$RIBOINSTALLDIR/vecscreen_plus_taxonomy" \
    PERL5LIB="$RIBOSCRIPTSDIR":"$RIBOSEQUIPDIR":"$VECPLUSDIR":"$PERL5LIB" \
    PATH="$RIBOSCRIPTSDIR":"$RRNASENSORDIR":"$PATH"

ENV \
    PATH="$RIBOSCRIPTSDIR":"$RRNASENSORDIR":"$PATH" \
    PATH="$RNA/traveler/bin:$RIBODIR:$RIBOEASELDIR:$RIBOINFERNALDIR:/rna/rscape/bin:/rna/jiffy-infernal-hmmer-scripts:/rna/r2dt:/rna/RNAstructure/exe:$PATH" \
    PERL5LIB="$RIBOSCRIPTSDIR:$RIBOSEQUIPDIR:$VECPLUSDIR:$BIOEASELDIR:$RIBODIR:$EPNOPTDIR:$EPNOFILEDIR:$EPNTESTDIR:$PERL5LIB:/usr/lib/tRNAscan-SE"

COPY --from=rscape /rna/rscape/ ${RNA}/rscape

COPY --from=trnascan-se /install/ /usr
RUN sed -i -e 's|/install/bin|/usr/bin|' -e 's|/install/lib|/usr/lib|' /usr/bin/tRNAscan-SE.conf

COPY --from=bio-easel /install/usr/* /usr
COPY --from=bio-easel /install/Bio-Easel $RNA/Bio-Easel
COPY --from=bio-easel /usr/local/share/perl/5.36.0/ /usr/share/perl/5.36/

COPY --from=scripts /install/ ${RNA}/jiffy-infernal-hmmer-scripts/

COPY --from=ribovore-infernal-easel /install/ ${RNA}/ribovore
COPY --from=rnartist /rnartist.jar /usr/local/bin/

COPY --from=traveler /install/traveler /usr/local/bin/traveler
COPY --from=traveler /install/utils/* ${RNA}/traveler/utils/
COPY --from=traveler "/usr/local/lib64/libstdc++.so.6.0.32" "/usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.30"
COPY --from=traveler "/usr/local/lib64/libstdc++.so.6.0.32" "/usr/lib/aarch64-linux-gnu/libstdc++.so.6.0.30"
