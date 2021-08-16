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

## Reference genome and other resources

The pipeline currently expects a number of files to exist on your system; we will consider including some functionality to build these on the flye if there is interest. Otherwise, the following applies:

1) all references are located under one common path on your shared file system, indicated by the variable `params.gatk_bundle_path` (this is a left-over from when we still used GATK). 

2) Under the common path, each assembly has one folder, being one or more of:

* hg19
* b37
* hg38
* hg38_no_alt

Each folder will contain a genome assembly in fasta format (check resources.config for details) and a list of dbSNP variants in VCF format. One way to get these two things is the GATK resource bundle - but note that `hg38_no_alt` is 
a custom creation, so is not included with the GATK bundle. To recreate it, simply download the ALT free GRCh38 assembly from NCBI and index it as described below. You can take the dbSNP file for hg38 from GATK, but need to re-header it with the sequence dictionary of your custom hg38 without ALT assembly.

For DeepVariant to work, the assembly has to be indexed in multiple ways:

fasta.fai `samtools faidx genome.fasta`
fasta.gz `bgzip -c -i genome.fasta`
fasta.gz.gzi
fasta.gz.fai `samtools faidx genome.fasta.gz`

Finally, the assembly-specific folder will contain a sub-folder called bwa2, containing (again) the (sym-linked) genome sequence and all the bwa2 relevant indices (created with `bwa index genome.fasta`).

Granted, the config is a bit clunky since it is more or less hard-coded for our systems, requiring you to make a number of custom adjustments on your end. If there is interest, we will consider making this a bit more flexible. 


