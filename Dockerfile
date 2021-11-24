FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for IKMB exome pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/exome-seq-3.1/bin:/opt/genesplicer/sources/:$PATH

RUN apt-get -y update && apt-get -y install make wget

RUN mkdir -p /opt && cd /opt && wget ftp://ftp.ccb.jhu.edu/pub/software/genesplicer/GeneSplicer.tar.gz && tar -xvf GeneSplicer.tar.gz && rm *.tar.gz && mv GeneSplicer genesplicer \
	&& cd genesplicer/sources/ && make 

