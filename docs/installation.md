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

This would add a new profile, called `your_profile` which uses (and expects) singularity to provide all software. 

`base.config` Basic settings about resource usage for the individual pipeline stages. 

`resources.config` Gives information about the files that are to be used during analysis for the individual human genome assemblies. 

`your_cluster.config` Specifies which sort of resource manager to use and where to find the resource bundle on your cluster file system (see below).

Alternatively, you can provide the config file during execution like so:

```
nextflow run ikmb/exome-seq -c my.config <other options here>
```

## Reference genome and other resources

The pipeline can build all relevant resources for this pipeline automatically. For a list of supported assemblies, see below. 

To do this, run:

```
nextflow run ikmb/exome-seq -c my.config --build_references --assembly GRCh38_no_alt --outdir /path/to/outdir
```

Note that `--outdir`specifies the root directory in which the reference files will be stored. You can provide this path to the main pipeline as `--genomes_base` or set this variable in your site specific config file. The pipeline should now be able to run with these pre-built references without requiring further downloads. 

Allowed reference assemblies are:

* GRCh38 (patch 1, with decoys and masked PAR regions - see [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) )
* GRCh38_no_alt (patch 1, no ALT contigs, with decoys and masked PAR regions - see [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) )
* GRCh38_g1k (patch1, as used by the 1000 genomes consortium and the Illumina Dragen system - see [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/))
* GRCh38_p14 (patch 14 without further modifications, see [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/) )
* GRCh38_no_alt_p14 (patch 14 without ALT contigs, see [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/) )
* hg38 (the BROAD version of GRCh38, part of the GATK bundle, see [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) )
* CHM13v2 (The telomere-to-telomere reference with the fully assembled Y chromosome, masked PAR regions and the mitochondrial genome rCRS, see [here](https://github.com/marbl/CHM13) )

If you plan on using effect prediction via VEP, you will also need to install the cache and plugins:

```
nextflow run ikmb/exome-seq -c my.config --build_references --assembly GRCh38_no_alt --outdir /path/to/outdir --tools vep
```

You only have to do this for one of the assemblies, as the VEP cache is the same for all supported assembly versions. 

### VEP reference files

Apart from the cache and plugin directory, VEP accepts additional reference files that can be used to annotate variants. These can be set in your local config file, although we cannot provide information on how to create them. Some support is available from the EnsEMBL VEP [plugin](https://github.com/Ensembl/VEP_plugins) code. 

- `dbnsfp_db`
- `dbscsnv_db`
- `cadd_snps`
- `cadd_indels`
- `vep_mastermind`


