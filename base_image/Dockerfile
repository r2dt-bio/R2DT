FROM gcc:13

ENV tRNAscan_VER=2.0.11 \
    BioEasel_COMMIT=e7fae0ab43fc47058183b71ff498ee0d0d2de6a7 \
    BioEasel_VER=0.09 \
# Note: Changing version from 2.4.x to other versions may require changes to the ViennaRNA installation code below
    ViennaRNA_VER=2.6.3 \
    autoconf_VER=2.71 \
    RNA=/rna


RUN apt-get update && apt-get install -q -y \
    gzip \
    less \
    moreutils \
    python3 \
    python3-pip \
    time \
    vim \
    wget \
    && rm -rf /var/lib/{apt,dpkg,cache,log}

RUN mkdir $RNA
WORKDIR $RNA


RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba && chmod +X bin/micromamba && mkdir -p /root/micromamba
RUN eval "$(./bin/micromamba shell hook -s posix)" && \
    micromamba activate && \
    micromamba config append channels defaults && \
    micromamba config append channels bioconda && \
    micromamba config append channels conda-forge && \
    micromamba config set channel_priority strict && \
    micromamba install viennarna


# Install R-scape
RUN wget http://rivaslab.org/software/rscape/rscape.tar.gz && \
    tar -xzf rscape.tar.gz && \
    mv rscape_* rscape && \
    rm rscape.tar.gz && \
    cd rscape && \
    cp /usr/share/automake-1.16/config.guess $RNA/rscape/lib/R2R/R2R-current/config.guess && \
    ./configure -q && make -j 4 && make install


# Install tRNAScan-SE
RUN \
    wget https://github.com/UCSC-LoweLab/tRNAscan-SE/archive/v${tRNAscan_VER}.tar.gz && \
    tar -xzf v${tRNAscan_VER}.tar.gz && \
    rm v${tRNAscan_VER}.tar.gz && \
    cd tRNAscan-SE-${tRNAscan_VER} && \
    ./configure -q && make -j 4 && make install

# Install Bio-Easel
RUN \
    cpan install Inline Inline::C && \
    git clone https://github.com/nawrockie/Bio-Easel.git && \
    cd Bio-Easel && \
    git checkout ${BioEasel_COMMIT} && \
    mkdir src && \
    cd src && \
    wget -O easel-Bio-Easel-${BioEasel_VER}.zip https://github.com/EddyRivasLab/easel/archive/Bio-Easel-${BioEasel_VER}.zip && \
    unzip easel-Bio-Easel-${BioEasel_VER}.zip && \
    mv easel-Bio-Easel-${BioEasel_VER} easel && \
    rm easel-Bio-Easel-${BioEasel_VER}.zip && \
    cd .. && \
    perl Makefile.PL; make -j 4; make install

# Install jiffy infernal hmmer scripts
RUN \
    git clone https://github.com/nawrockie/jiffy-infernal-hmmer-scripts.git && \
    cd jiffy-infernal-hmmer-scripts && \
    git checkout 31d6c3b826d432a30620507830749cee58e15e68 && \
    echo '#!/usr/bin/env perl' | cat - ali-pfam-sindi2dot-bracket.pl | sponge ali-pfam-sindi2dot-bracket.pl && \
    chmod +x $RNA/jiffy-infernal-hmmer-scripts/*.pl

## Install ViennaRNA
#RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_6_x/ViennaRNA-${ViennaRNA_VER}.tar.gz && \
#    tar -zxvf ViennaRNA-${ViennaRNA_VER}.tar.gz && \
#    cd ViennaRNA-${ViennaRNA_VER} && \
#    ./configure -q --with-python3 && \
#    make -j 4 && \
#    make install && \
#    cd $RNA && \
#    rm ViennaRNA-${ViennaRNA_VER}.tar.gz



# Install develop branch of easel
RUN git clone https://github.com/EddyRivasLab/easel && \
    cd easel && \
    autoconf && \
    ./configure -q && \
    make -j 4

# Install Traveler
RUN git clone https://github.com/cusbg/traveler && \
    cd traveler && \
    git checkout v3.0.0 & \
    cd src && make -j 4 build

## Install autoconf required by Infernal
#RUN \
#    wget http://ftp.gnu.org/gnu/autoconf/autoconf-${autoconf_VER}.tar.gz && \
#    tar xzvf autoconf-${autoconf_VER}.tar.gz && \
#    cd autoconf-${autoconf_VER} && \
#    ./configure -q && \
#    make -j 4 && make install && \
#    cd $RNA && \
#    rm autoconf-${autoconf_VER}.tar.gz && \
#    rm -Rf autoconf-${autoconf_VER}


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

# Setup environmental variables
ENV \
    RIBOINSTALLDIR="$RNA/ribovore" \
    RIBOEASELDIR="$RNA/ribovore/infernal/easel/miniapps" \
    RIBOTIMEDIR=/usr/bin \
    RIBODIR="$RNA/ribovore/ribovore" \
    BIOEASELDIR="$RNA/Bio-Easel/blib/lib:$RNA/Bio-Easel/blib/arch:$RNA/Bio-Easel:$RNA/Bio-Easel/lib" \
    LC_ALL="C.UTF-8" LANG="C.UTF-8" \
    DATAPATH="/rna/RNAstructure/data_tables/" \
    PYTHONPATH="$PYTHONPATH:/rna/ViennaRNA-${ViennaRNA_VER}/interfaces/Python3"
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

WORKDIR ${RNA}/r2dt

ENTRYPOINT ["/bin/bash"]
