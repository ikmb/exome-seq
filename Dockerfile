FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for IKMB exome pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/exome-seq-1.3/bin:$PATH
RUN mkdir -p /opt/vep/ && cd /opt/vep/ && git clone git@github.com:Ensembl/VEP_plugins.git plugins && cd plugins &&\
git checkout release/99
