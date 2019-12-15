FROM r-base:3.6.1

RUN mkdir /scripts /soft

COPY /scripts/. /scripts
COPY /soft/. /soft

# Update sources list, install wget, unzip, python, python3
RUN apt-get update && apt-get install --yes wget unzip python python3

# Working directory in /bin
WORKDIR /bin

# Download conda
RUN wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O /bin/miniconda.sh
RUN chmod 0755 /bin/miniconda.sh
RUN /bin/miniconda.sh -b -p /bin/conda
ENV PATH="/bin/conda/bin:$PATH"
RUN rm miniconda.sh

# Update conda
RUN conda update conda --yes
RUN conda install conda-build --yes

# Download TopHat2
RUN wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.9.Linux_x86_64.tar.gz
RUN tar zxvf tophat-2.0.9.Linux_x86_64.tar.gz
RUN rm tophat-2.0.9.Linux_x86_64.tar.gz
ENV PATH=$PATH:/bin/tophat-2.0.9.Linux_x86_64

# Download Bowtie2
RUN wget --default-page=bowtie2-2.3.4.1-linux-x86_64.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip
RUN unzip bowtie2-2.1.0-linux-x86_64.zip
RUN rm bowtie2-2.1.0-linux-x86_64.zip
ENV PATH=$PATH:/bin/bowtie2-2.1.0

# Download Bowtie
RUN wget --default-page=bowtie-1.0.0-linux-x86_64.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.0.0/bowtie-1.0.0-linux-x86_64.zip
RUN unzip bowtie-1.0.0-linux-x86_64.zip
RUN rm bowtie-1.0.0-linux-x86_64.zip
ENV PATH=$PATH:/bin/bowtie-1.0.0

# Download miRDeep2
RUN conda install -c bioconda mirdeep2=2.0.1.2 --yes
ENV PATH=$PATH:/bin/conda/pkgs/mirdeep2-2.0.1.2-0/bin

# Download CIRCExplorer2
RUN conda install -c bioconda circexplorer2=2.3.6 --yes
ENV PATH=$PATH:/bin/conda/pkgs/circexplorer2-2.3.6-py_0/site-packages/

# Downloand miRanda
RUN wget http://cbio.mskcc.org/microrna_data/miRanda-aug2010.tar.gz
RUN tar -xvzf miRanda-aug2010.tar.gz
WORKDIR /bin/miRanda-3.3a/
RUN ./configure
RUN make install
WORKDIR /bin/

# Download samtools
RUN conda install -c bioconda samtools=0.1.19

RUN Rscript /soft/install_packages.R
RUN Rscript /soft/test_installation.R






