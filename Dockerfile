FROM gcc:6

RUN apt-get update && apt-get install -y moreutils python3 python3-pip less vim

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
    wget http://trna.ucsc.edu/software/trnascan-se-2.0.5.tar.gz && \
    tar -xvzf trnascan-se-2.0.5.tar.gz && \
    rm trnascan-se-2.0.5.tar.gz && \
    cd tRNAscan-SE-2.0 && \
    ./configure && make && make install
# Make sure tRNAScan-SE can find Infernal
RUN \
    ln -s /rna/infernal-1.1.2/src/cmsearch /usr/local/bin/cmsearch && \
    ln -s /rna/infernal-1.1.2/src/cmscan /usr/local/bin/cmscan

# Install jiffy infernal hmmer scripts
RUN \
    git clone https://github.com/nawrockie/jiffy-infernal-hmmer-scripts.git && \
    cd jiffy-infernal-hmmer-scripts && \
    git checkout 45d4937385a6b694eac2d7d538e131b59527ce06
RUN \
    cd jiffy-infernal-hmmer-scripts && \
    echo '#!/usr/bin/env perl' | cat - ali-pfam-sindi2dot-bracket.pl | sponge ali-pfam-sindi2dot-bracket.pl
RUN chmod +x $RNA/jiffy-infernal-hmmer-scripts/ali-pfam-sindi2dot-bracket.pl

# Install ribovore
RUN git clone https://github.com/nawrockie/epn-ofile.git && cd epn-ofile && git fetch && git fetch --tags && git checkout ribovore-0.38
RUN git clone https://github.com/nawrockie/epn-options.git && cd epn-options && git fetch && git fetch --tags && git checkout ribovore-0.38
RUN git clone https://github.com/nawrockie/epn-test.git && cd epn-test && git fetch && git fetch --tags && git checkout ribovore-0.38
RUN git clone https://github.com/nawrockie/ribovore.git && cd ribovore && git checkout auto-traveler

# Install Traveler
RUN git clone https://github.com/davidhoksza/traveler.git && cd traveler && git checkout a87234b179ebea6cac213ffd9e675d580dd60885
RUN cd $RNA/traveler/src && make build

COPY examples examples/

# Install python dependencies
ADD requirements.txt $RNA/auto-traveler/requirements.txt
RUN pip3 install -r $RNA/auto-traveler/requirements.txt

# Setup environmental variables
ENV RIBODIR="$RNA/ribovore" RIBOINFERNALDIR="$RNA/infernal-1.1.2/bin" RIBOEASELDIR="$RNA/infernal-1.1.2/bin"
ENV EPNOPTDIR="$RNA/epn-options" EPNOFILEDIR="$RNA/epn-ofile" EPNTESTDIR="$RNA/epn-test"
RUN apt-get update && apt-get install time
ENV RIBOTIMEDIR="/usr/bin"
ENV PERL5LIB="$RIBODIR:$EPNOPTDIR:$EPNOFILEDIR:$EPNTESTDIR:$PERL5LIB"
ENV LC_ALL="C.UTF-8" LANG="C.UTF-8"
ENV PATH="$RNA/traveler/bin:$RIBODIR:$RIBOINFERNALDIR:$PATH"
ENV PATH="/rna/rscape/bin:$PATH"
ENV PATH="/rna/jiffy-infernal-hmmer-scripts/:$PATH"
ENV PATH="/rna/RNAstructure/exe:$PATH" DATAPATH="/rna/RNAstructure/data_tables/"
ENV PATH="/rna/auto-traveler:$PATH"

ENTRYPOINT ["/bin/bash"]
