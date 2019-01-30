FROM gcc:4.9

RUN apt-get update
RUN apt-get install -y moreutils

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

# Install ribotyper-v1
RUN git clone https://github.com/nawrockie/epn-ofile.git && cd epn-ofile && git checkout c34244b2b9e0719c45d964cc08c147aa353532e8
RUN git clone https://github.com/nawrockie/epn-options.git && cd epn-options && git checkout 7acc13384aedbd5efee9a62fcde71d075072b6a6
RUN git clone https://github.com/nawrockie/epn-test.git && cd epn-test && git checkout f4a8a60153906e61bc458fa734ec7070eadf76f9
RUN git clone https://github.com/nawrockie/ribotyper-v1.git && cd ribotyper-v1 && git checkout 4cd7fe30f402edfa4669383a46d603c60ba6f608

# Install Traveler
RUN git clone https://github.com/davidhoksza/traveler.git && cd traveler && git checkout 0912ed5daab09bb3c38630efaf3643ea38b02dbe
RUN cd $RNA/traveler/src && make build

# Setup environmental variables
ENV RIBODIR="$RNA/ribotyper-v1" RIBOINFERNALDIR="$RNA/infernal-1.1.2/bin" RIBOEASELDIR="$RNA/infernal-1.1.2/bin"
ENV EPNOPTDIR="$RNA/epn-options" EPNOFILEDIR="$RNA/epn-ofile" EPNTESTDIR="$RNA/epn-test"
ENV PERL5LIB="$RIBODIR:$EPNOPTDIR:$EPNOFILEDIR:$EPNTESTDIR:$PERL5LIB"
ENV PATH="$RNA/traveler/bin:$RIBODIR:$RIBOINFERNALDIR:$PATH"

# Install jiffy infernal hmmer scripts
RUN \
    git clone https://github.com/nawrockie/jiffy-infernal-hmmer-scripts.git && \
    cd jiffy-infernal-hmmer-scripts && \
    git checkout 45d4937385a6b694eac2d7d538e131b59527ce06
RUN \
    cd jiffy-infernal-hmmer-scripts && \
    echo '#!/usr/bin/env perl' | cat - ali-pfam-sindi2dot-bracket.pl | sponge ali-pfam-sindi2dot-bracket.pl
RUN chmod +x $RNA/jiffy-infernal-hmmer-scripts/ali-pfam-sindi2dot-bracket.pl

COPY examples examples/

# Install RNAStructure
RUN \
    wget http://rna.urmc.rochester.edu/Releases/current/RNAstructureSource.tgz && \
    tar -xvzf RNAstructureSource.tgz && \
    rm RNAstructureSource.tgz && \
    cd RNAstructure && \
    make all

# Install pip
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python get-pip.py

# Install python dependencies
ADD requirements.txt $RNA/auto-traveler/requirements.txt
RUN pip install -r $RNA/auto-traveler/requirements.txt

ENV PATH="/rna/jiffy-infernal-hmmer-scripts/:$PATH"
ENV PATH="/rna/RNAstructure/exe:$PATH" DATAPATH="/rna/RNAstructure/data_tables/"

ENTRYPOINT ["/bin/bash"]
