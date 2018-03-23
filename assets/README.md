![](images/ikmb_bfx_logo.png)

# IKMB Diagnostic Exome Pipeline - Singularity Container

Singularity (http://singularity.lbl.gov/) is a container-framework similar to Docker. In order to build the container, you will  need the following:

- A Linux server
- Root access
- The bootstrap file included with this code
- A working installation of singularity (tested with 2.4)

Depending on the flavor of your host linux, you may also need to install additional dependencies - you will be informed of this when you try building the container. 

## Building the container

`sudo singularity image.create -s 4200 container.img`

This builds an empty container with a fixed size of 4.2GB.

`sudo singularity build container.img ikmb-diagnostic-exomes.definition`

Et voila. 

As a final step, you will need to create a config file (see config/template for a blank template). Specifically, the config file requires information about where
the container is located as well as any and all settings related to your resource management system (suxh as Slurm, LSF etc). If you are unsure, you can check the other config files, such as rzcluster.config which is configured for a Slurm cluster consisting of nodes with 16 cores and 128GB Ram. 

## Reference data

This pipeline uses a number of reference data sets, including genome assemblies and data indices. The location of these data are specified in the cluster config file (see config/rzcluster.config for an example). The data used by us comes from the BROAD resource bundle; simply copy it to some location on your cluster and set the relevant paths in your config file.
