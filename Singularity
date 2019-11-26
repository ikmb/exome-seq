Bootstrap:docker
From:continuumio/anaconda

%labels
    MAINTAINER Marc Hoeppner <m.hoeppner@ikmb.uni-kiel.de>
    DESCRIPTION Singularity image containing all requirements for the exome-seq pipeline
    VERSION 1.2

%environment
    PATH=/opt/conda/envs/exome-seq-1.2/bin:$PATH
    export PATH

%files
    environment.yml /

%post

    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a

    mkdir -p /ifs
    apt-get -y install procps

