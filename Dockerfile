FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for IKMB exome pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/exome-seq-1.1/bin:$PATH

RUN apt-get -y update && apt-get -y install texlive texlive-latex-extra texlive-generic-extra texlive-xetex texlive-latex-recommended texlive-luatex
