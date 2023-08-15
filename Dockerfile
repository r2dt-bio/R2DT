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

# Install Infernal
RUN \
    cd $RNA && \
    curl -OL http://eddylab.org/infernal/infernal-1.1.2.tar.gz && \
    tar -xvzf infernal-1.1.2.tar.gz && \
    cd infernal-1.1.2 && \
    ./configure --prefix=$RNA/infernal-1.1.2 && \
    make && \
    make install && \
    cd easel && \
    make install && \
    cd $RNA && \
    rm infernal-1.1.2.tar.gz

# Install R-scape
RUN wget http://eddylab.org/software/rscape/rscape.tar.gz && \
    tar -xvzf rscape.tar.gz && \
    mv rscape_* rscape && \
    rm rscape.tar.gz && \
    cd rscape && \
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
# Make sure tRNAScan-SE can find Infernal
RUN \
    ln -s /rna/infernal-1.1.2/src/cmsearch /usr/local/bin/cmsearch && \
    ln -s /rna/infernal-1.1.2/src/cmscan /usr/local/bin/cmscan

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

# Install ribovore
RUN git clone https://github.com/nawrockie/epn-ofile.git && cd epn-ofile && git fetch && git fetch --tags && git checkout ribovore-0.40
RUN git clone https://github.com/nawrockie/epn-options.git && cd epn-options && git fetch && git fetch --tags && git checkout ribovore-0.40
RUN git clone https://github.com/nawrockie/epn-test.git && cd epn-test && git fetch && git fetch --tags && git checkout ribovore-0.40
RUN git clone https://github.com/ncbi/ribovore.git && cd ribovore && git fetch && git fetch --tags && git checkout ribovore-0.40

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

# Install python dependencies
ADD . /rna/r2dt
ADD requirements.txt $RNA/r2dt/requirements.txt
RUN pip3 install -r $RNA/r2dt/requirements.txt

# Setup environmental variables
ENV RIBODIR="$RNA/ribovore" \
    RIBOINFERNALDIR="$RNA/infernal-1.1.2/bin" \
    RIBOEASELDIR="$RNA/infernal-1.1.2/bin" \
    EPNOPTDIR="$RNA/epn-options" \
    EPNOFILEDIR="$RNA/epn-ofile" \
    EPNTESTDIR="$RNA/epn-test" \
    RIBOTIMEDIR="/usr/bin" \
    BIOEASELDIR="$RNA/Bio-Easel/blib/lib:$RNA/Bio-Easel/blib/arch:$RNA/Bio-Easel:$RNA/Bio-Easel/lib" \
    LC_ALL="C.UTF-8" LANG="C.UTF-8" \
    DATAPATH="/rna/RNAstructure/data_tables/" \
    PYTHONPATH="$PYTHONPATH:/rna/ViennaRNA-2.4.18/interfaces/Python3"
ENV PATH="$RNA/traveler/bin:$RIBODIR:$RIBOINFERNALDIR:/rna/rscape/bin:/rna/jiffy-infernal-hmmer-scripts:/rna/r2dt:/rna/RNAstructure/exe:$PATH" \
    PERL5LIB="$BIOEASELDIR:$RIBODIR:$EPNOPTDIR:$EPNOFILEDIR:$EPNTESTDIR:$PERL5LIB"
WORKDIR /rna/r2dt

ENTRYPOINT ["/bin/bash"]
