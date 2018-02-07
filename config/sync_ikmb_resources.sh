#!/bin/bash

# Remote information
USER="sukmb352"
CONNECTION="${USER}@ikmbhead.rz.uni-kiel.de"

# Local information
STORAGE="/storage2"
REFERENCES="${STORAGE}/references"
DATABASES="${STORAGE}/databases"
IMAGES="${STORAGE}/images"

# Staging
mkdir -p $IMAGES
mkdir -p $DATABASES
mkdir -p $DATABASES/annovar
mkdir -p $REFERENCES
mkdir -p $REFERENCES/exomes

# Data transfer
rsync --exclude 'bundle/hg38' --exclude 'bundle/test' --exclude 'bundle/2.8/b36' --exclude 'bundle/2.8/hg18' --exclude 'bundle/2.8/exampleFASTA' --exclude 'bundle/2.8/minimalhg19' -av $CONNECTION:/ifs/data/nfs_share/ikmb_repository/references/gatk $REFERENCES
rsync -av --exclude '*.bam*' $CONNECTION:/ifs/data/nfs_share/ikmb_repository/references/exomes/calibration_exomes_resorted $REFERENCES/exomes
rsync -av $CONNECTION:/ifs/data/nfs_share/ikmb_repository/references/exomes/nextera_exome_target_2017 $REFERENCES/exomes
rsync -av $CONNECTION:/ifs/data/nfs_share/ikmb_repository/references/exomes/idt_xgen $REFERENCES/exomes

rsync -av $CONNECTION:/ifs/data/nfs_share/ikmb_repository/software/annovar/2017-09-20/humandb $DATABASES/annovar

rsync -av $CONNECTION:/ifs/data/nfs_share/ikmb_repository/references/singularity/images/ikmb-diagnostic-exomes.img $IMAGES
