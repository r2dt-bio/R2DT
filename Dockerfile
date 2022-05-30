FROM gcc:10

RUN apt-get update && apt-get install -y moreutils python3 python3-pip gzip less wget time vim

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

# Install RNAStructure
RUN \
    wget http://rna.urmc.rochester.edu/Releases/current/RNAstructureSource.tgz && \
    tar -xvzf RNAstructureSource.tgz && \
    rm RNAstructureSource.tgz && \
    cd RNAstructure && \
    make all

# Install tRNAScan-SE
RUN \
    wget https://github.com/UCSC-LoweLab/tRNAscan-SE/archive/v2.0.5.tar.gz && \
    tar -xvzf v2.0.5.tar.gz && \
    rm v2.0.5.tar.gz && \
    cd tRNAscan-SE-2.0.5 && \
    ./configure && make && make install
# Make sure tRNAScan-SE can find Infernal
RUN \
    ln -s /rna/infernal-1.1.2/src/cmsearch /usr/local/bin/cmsearch && \
    ln -s /rna/infernal-1.1.2/src/cmscan /usr/local/bin/cmscan

# Install Bio-Easel
RUN \
    cpan install Inline && \
    cpan install Inline::C
RUN \
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
    perl Makefile.PL; make; make test; make install

# Install Python 3.6
RUN \
    mkdir python36 && \
    wget https://www.python.org/ftp/python/3.6.11/Python-3.6.11.tgz && \
    tar -xvf Python-3.6.11.tgz && \
    cd Python-3.6.11 && \
    ./configure --prefix=$RNA/python36/ && \
    make && make install && \
    cd .. && \
    rm -Rf Python-3.6.11

# Install jiffy infernal hmmer scripts
RUN \
    git clone https://github.com/nawrockie/jiffy-infernal-hmmer-scripts.git && \
    cd jiffy-infernal-hmmer-scripts && \
    git checkout 23b78b30b49b5255bae2cebb3f96fd3a147059a6
RUN \
    cd jiffy-infernal-hmmer-scripts && \
    echo '#!/usr/bin/env perl' | cat - ali-pfam-sindi2dot-bracket.pl | sponge ali-pfam-sindi2dot-bracket.pl && \
    chmod +x $RNA/jiffy-infernal-hmmer-scripts/*.pl

# Install ribovore
RUN git clone https://github.com/nawrockie/epn-ofile.git && cd epn-ofile && git fetch && git fetch --tags && git checkout ribovore-0.40
RUN git clone https://github.com/nawrockie/epn-options.git && cd epn-options && git fetch && git fetch --tags && git checkout ribovore-0.40
RUN git clone https://github.com/nawrockie/epn-test.git && cd epn-test && git fetch && git fetch --tags && git checkout ribovore-0.40
RUN git clone https://github.com/ncbi/ribovore.git && cd ribovore && git fetch && git fetch --tags && git checkout ribovore-0.40

#Install ViennaRNA
RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.18.tar.gz && \
    tar -zxvf ViennaRNA-2.4.18.tar.gz && \
    cd ViennaRNA-2.4.18 && \
    ./configure --with-python3 && \
    make && \
    make install

# Install Traveler
RUN git clone https://github.com/cusbg/traveler && cd traveler && git checkout 5e363fe3079ce5b80103889b9c9e213e1d2a16ff && cd src && make build

COPY examples examples/

# Install python dependencies
ADD . /rna/r2dt
ADD requirements.txt $RNA/r2dt/requirements.txt
RUN pip3 install -r $RNA/r2dt/requirements.txt

# Setup environmental variables
ENV RIBODIR="$RNA/ribovore" RIBOINFERNALDIR="$RNA/infernal-1.1.2/bin" RIBOEASELDIR="$RNA/infernal-1.1.2/bin"
ENV EPNOPTDIR="$RNA/epn-options" EPNOFILEDIR="$RNA/epn-ofile" EPNTESTDIR="$RNA/epn-test"
ENV RIBOTIMEDIR="/usr/bin"
ENV BIOEASELDIR="$RNA/Bio-Easel/blib/lib:$RNA/Bio-Easel/blib/arch:$RNA/Bio-Easel:$RNA/Bio-Easel/lib"
ENV PERL5LIB="$BIOEASELDIR:$RIBODIR:$EPNOPTDIR:$EPNOFILEDIR:$EPNTESTDIR:$PERL5LIB"
ENV LC_ALL="C.UTF-8" LANG="C.UTF-8"
ENV PATH="$RNA/traveler/bin:$RIBODIR:$RIBOINFERNALDIR:$PATH"
ENV PATH="/rna/rscape/bin:$PATH"
ENV PATH="/rna/jiffy-infernal-hmmer-scripts:$PATH"
ENV PATH="/rna/RNAstructure/exe:$PATH" DATAPATH="/rna/RNAstructure/data_tables/"
ENV PATH="/rna/r2dt:$PATH"
ENV PYTHONPATH="$PYTHONPATH:/rna/ViennaRNA-2.4.18/interfaces/Python3"
WORKDIR /rna/r2dt

ENTRYPOINT ["/bin/bash"]
