FROM jupyter/scipy-notebook

USER root

# Install CMake, GSL and xcas (for giac)

RUN rm -rf /var/lib/apt/lists/* && apt-get clean
RUN  apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential cmake libgsl-dev xcas 

# Install SLiM
RUN wget http://benhaller.com/slim/SLiM.zip 
RUN unzip SLiM.zip
RUN mkdir SLiM/build 
RUN cd SLiM/build && cmake .. && make install -j 4
RUN rm -fR SLiM*

RUN conda install --quiet --yes \
    zarr scikit-allel \
    more-itertools tqdm sympy networkx psutil pandas \
    docopt pytables tabulate htop \
    && conda clean -tipsy 
 
RUN conda install --quiet --yes \
    -c bioconda pysam \
    && conda clean -tipsy 

RUN pip install --upgrade --pre tskit==0.2.0a3
RUN pip install msprime tsinfer pyslim tszip

# nbgitpuller to pull in data files, images etc.
RUN pip install nbgitpuller
RUN jupyter serverextension enable --py nbgitpuller --sys-prefix

RUN cd /tmp && wget https://www.dropbox.com/s/hhz959alniylfpm/hmel.chr18.tar.gz \
    && mkdir -p /data/hmel.chr18 && cd /data/hmel.chr18 && tar -zxvf /tmp/hmel.chr18.tar.gz && \
    rm -fR /tmp/hmel.chr18.tar.gz

RUN fix-permissions /home/jovyan
USER $NB_UID
