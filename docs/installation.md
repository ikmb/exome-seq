# Installation

## At the IKMB

If you are at the IKMB, you will not have to do anything to make this run, it is all pre-configured for our compute system(s).
 
## Site-specific config file

This pipeline requires a site-specific configuration file to be able to talk to your local cluster or compute infrastructure. Nextflow supports a wide
range of such infrastructures, including Slurm, LSF and SGE - but also Kubernetes and AWS. For more information, see [here](https://www.nextflow.io/docs/latest/executor.html).

Please see conf/medcluster.config for an example of how to configure this pipeline for a Slurm queue.

All software is provided through Docker containers - this requires for your compute system to run either Docker or Singularity (more common on HPC systems). Details on how to specify singularity as your container engine are provided in the config file for our medcluster (medcluster.config).

With this information in place, you will next have to create an new site-specific profile for your local environment in `nextflow.config` using the following format:

```

profiles {
	
	your_profile {
		includeConfig 'conf/base.config'
		includeConfig 'conf/your_cluster.config'
		includeConfig 'conf/resources.config'
	}
}

```

This would add a new profile, called `your_profile` which uses (and expects) conda to provide all software. 

`base.config` Basic settings about resource usage for the individual pipeline stages. 

`resources.config` Gives information about the files that are to be used during analysis for the individual human genome assemblies. 

`your_cluster.config` Specifies which sort of resource manager to use and where to find the GATK resource bundle on your cluster file system (see below).

## GATK resource bundle

This pipeline uses GATK for variant calling and as such it requires a fairly large number of specific reference files, including a genome assembly and specifically matched variant references for annotation and filtering as well as exome kit definitions. 

The by far easiest way to ensure that all these files are available is to download the GATK [resource bundle](https://software.broadinstitute.org/gatk/download/bundle) to your local cluster. 
Nextflow also understands S3 buckets as file locations, if that fits your use case better. 

Once downloaded, you will have to decompress the genome assembly and ensure that the required BWA index files are present (they should be included). The root directory, i.e. the one directory containing the various genome assembly sub folders, must then be added into your site-specific config file `your_cluster.config` using:

`gatk_bundle_path = "/path/to/resource_bundle/" `

