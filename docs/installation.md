# Installation

## Site-specific config file

This pipeline requires a site-specific configuration file to be able to talk to your local cluster or compute infrastructure. Nextflow supports a wide
range of such infrastructures, including Slurm, LSF and SGE - but also Kubernetes and AWS. For more information, see [here](https://www.nextflow.io/docs/latest/executor.html).

Please see conf/rzcluster.config for an example of how to configure this pipeline for a Slurm queue.

In addition to specifying a cluster environment, you will also have to decide how you wish the pipeline to provision the necessary software.

Option 1: Install all relevant tools via conda during start-up (this requires [conda](https://anaconda.org/)). 

Option 2: Pull a pre-configured [Docker](https://cloud.docker.com/u/ikmb/repository/docker/ikmb/exome-seq) container. 


With this information in place, you will next have to create an new site-specific profile for your local environment in `nextflow.config` using the following format:

```

profiles {
	
	your_profile {
		includeConfig 'conf/base.config'
		includeConfig 'conf/your_cluster.config'
		includeConfig 'conf/resources.config'
		includeConfig 'conf/conda.config'
	}
}

```

This would add a new profile, called `your_profile` which uses (and expects) conda to provide all software. 

`base.config` Basic settings about resource usage for the individual pipeline stages. 

`resources.config` Gives information about the files that are to be used during analysis for the individual human genome assemblies. 

`conda.config` Specifies how conda is to be used for software provisioning. 

`singularity.config` Specifies how to use Singularity for software provisioning.

`your_cluster.config` Specifies which sort of resource manager to use and where to find the GATK resource bundle on your cluster file system (see below).

## GATK resource bundle

This pipeline uses GATK for variant calling and as such it requires a fairly large number of specific reference files, including a genome assembly and specifically matched variant references for annotation and filtering as well as exome kit definitions. 

The by far easiest way to ensure that all these files are available is to download the GATK [resource bundle](https://software.broadinstitute.org/gatk/download/bundle) to your local cluster. 
Nextflow also understands S3 buckets as file locations, if that fits your use case better. 

Once downloaded, you will have to decompress the genome assembly and ensure that the required BWA index files are present (they should be included). The root directory, i.e. the one directory containing the various genome assembly sub folders, must then be added into your site-specific config file `your_cluster.config` using:

`gatk_bundle_path = "/path/to/resource_bundle/" `

