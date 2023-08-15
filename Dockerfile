FROM gcc:11

RUN apt-get update && apt-get install -y \
    gzip \
    less \
    moreutils \
    python3 \
    python3-pip \
    time \
    vim \
    wget \
    && rm -rf /var/lib/apt/lists/*

ENV RNA /rna

RUN mkdir $RNA

WORKDIR $RNA

# Install R-scape
RUN wget http://rivaslab.org/software/rscape/rscape.tar.gz && \
    tar -xvzf rscape.tar.gz && \
    mv rscape_* rscape && \
    rm rscape.tar.gz && \
    cd rscape && \
    cp /usr/share/automake-1.16/config.guess $RNA/rscape/lib/R2R/R2R-current/config.guess && \
    ./configure && make && make install

# Install RNAStructure - only needed for updating CRW templates
# RUN \
#     wget http://rna.urmc.rochester.edu/Releases/current/RNAstructureSource.tgz && \
#     tar -xvzf RNAstructureSource.tgz && \
#     rm RNAstructureSource.tgz && \
#     cd RNAstructure && \
#     make all

# Install tRNAScan-SE
RUN \
    wget https://github.com/UCSC-LoweLab/tRNAscan-SE/archive/v2.0.11.tar.gz && \
    tar -xvzf v2.0.11.tar.gz && \
    rm v2.0.11.tar.gz && \
    cd tRNAscan-SE-2.0.11 && \
    ./configure && make && make install

# Install Bio-Easel
RUN \
    cpan install Inline && \
    cpan install Inline::C && \
    git clone https://github.com/nawrockie/Bio-Easel.git && \
    cd Bio-Easel && \
    git checkout e7fae0ab43fc47058183b71ff498ee0d0d2de6a7 && \
    mkdir src && \
    cd src && \
    curl -k -L -o easel-Bio-Easel-0.09.zip https://github.com/EddyRivasLab/easel/archive/Bio-Easel-0.09.zip && \
    unzip easel-Bio-Easel-0.09.zip && \
    mv easel-Bio-Easel-0.09 easel && \
    rm easel-Bio-Easel-0.09.zip && \
    cd .. && \
    perl Makefile.PL; make; make install

# Install jiffy infernal hmmer scripts
RUN \
    git clone https://github.com/nawrockie/jiffy-infernal-hmmer-scripts.git && \
    cd jiffy-infernal-hmmer-scripts && \
    git checkout 31d6c3b826d432a30620507830749cee58e15e68 && \
    echo '#!/usr/bin/env perl' | cat - ali-pfam-sindi2dot-bracket.pl | sponge ali-pfam-sindi2dot-bracket.pl && \
    chmod +x $RNA/jiffy-infernal-hmmer-scripts/*.pl

# Install ViennaRNA
RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.18.tar.gz && \
    tar -zxvf ViennaRNA-2.4.18.tar.gz && \
    cd ViennaRNA-2.4.18 && \
    ./configure --with-python3 && \
    make && \
    make install && \
    cd $RNA && \
    rm ViennaRNA-2.4.18.tar.gz

# Install develop branch of easel
RUN git clone https://github.com/EddyRivasLab/easel && \
    cd easel && \
    git checkout develop && \
    autoconf && \
    ./configure && \
    make && \
    make check

# Install Traveler
RUN git clone https://github.com/cusbg/traveler && \
    cd traveler && \
    git checkout 0ee67cbb9c5aaa2f98340065fd047c9a8feea53e && \
    cd src && \
    make build

# Install autoconf required by Infernal
RUN \
    wget http://ftp.gnu.org/gnu/autoconf/autoconf-2.71.tar.gz && \
    tar xzvf autoconf-2.71.tar.gz && \
    cd autoconf-2.71 && \
    ./configure && \
    make && \
    make install && \
    cd $RNA && \
    rm autoconf-2.71.tar.gz && \
    rm -Rf autoconf-2.71

# Install Ribovore and Infernal
RUN \
    mkdir -p $RNA/ribovore && \
    cd $RNA/ribovore && \
    wget https://raw.githubusercontent.com/ncbi/ribovore/r2dt/infdev-install.sh && \
    chmod +x infdev-install.sh && \
    ./infdev-install.sh "linux"

# Make sure tRNAScan-SE can find Infernal
RUN \
    ln -s /rna/ribovore/infernal/src/cmsearch /usr/local/bin/cmsearch && \
    ln -s /rna/ribovore/infernal/src/cmscan /usr/local/bin/cmscan

# Install python dependencies
ADD . /rna/r2dt
ADD requirements.txt $RNA/r2dt/requirements.txt
RUN pip3 install -r $RNA/r2dt/requirements.txt

# Setup environmental variables
ENV \
    RIBOINSTALLDIR="$RNA/ribovore" \
    RIBOEASELDIR="$RNA/ribovore/infernal/easel/miniapps" \
    RIBOTIMEDIR=/usr/bin \
    RIBODIR="$RNA/ribovore/ribovore" \
    BIOEASELDIR="$RNA/Bio-Easel/blib/lib:$RNA/Bio-Easel/blib/arch:$RNA/Bio-Easel:$RNA/Bio-Easel/lib" \
    LC_ALL="C.UTF-8" LANG="C.UTF-8" \
    DATAPATH="/rna/RNAstructure/data_tables/" \
    PYTHONPATH="$PYTHONPATH:/rna/ViennaRNA-2.4.18/interfaces/Python3"
ENV \
    RIBOSCRIPTSDIR="$RIBOINSTALLDIR/ribovore" \
    RIBOINFERNALDIR="$RIBOINSTALLDIR/infernal/src" \
    RIBOSEQUIPDIR="$RIBOINSTALLDIR/sequip" \
    RIBOBLASTDIR="$RIBOINSTALLDIR/ncbi-blast/bin" \
    VECPLUSDIR="$RIBOINSTALLDIR/vecscreen_plus_taxonomy" \
    BLASTDB="$VECPLUSDIR/univec-files":"$RRNASENSORDIR":"$BLASTDB" \
    RRNASENSORDIR="$RIBOINSTALLDIR/rRNA_sensor"
ENV \
    PATH="$RIBOSCRIPTSDIR":"$RIBOBLASTDIR":"$RRNASENSORDIR":"$PATH" \
    PATH="$RNA/traveler/bin:$RIBODIR:$RIBOEASELDIR:$RIBOINFERNALDIR:/rna/rscape/bin:/rna/jiffy-infernal-hmmer-scripts:/rna/r2dt:/rna/RNAstructure/exe:$PATH" \
    PERL5LIB="$RIBOSCRIPTSDIR:$RIBOSEQUIPDIR:$VECPLUSDIR:$BIOEASELDIR:$RIBODIR:$EPNOPTDIR:$EPNOFILEDIR:$EPNTESTDIR:$PERL5LIB"
WORKDIR /rna/r2dt

ENTRYPOINT ["/bin/bash"]
